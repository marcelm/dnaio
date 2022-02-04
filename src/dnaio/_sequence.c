#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef
#include "_sequence.h"

static void 
SequenceRecord_dealloc(SequenceRecord *self) {
    Py_CLEAR(self->name);
    Py_CLEAR(self->sequence);
    Py_CLEAR(self->qualities);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

PyDoc_STRVAR(SequenceRecord__init____doc__,
"SequenceRecord(name, sequence, qualities = None)\n"
"--\n"
"\n"
"A sequencing read with read name/id and (optional) qualities\n"
"\n"
"For FASTA, the qualities attribute is None.\n" 
"For FASTQ, qualities is a string and it contains the qualities\n"
"encoded as ASCII(qual+33)\n."
"\n"
"  name\n"
"    The read description\n"
"  sequence\n"
"  qualities\n"
"\n"
);

static int
SequenceRecord__init__(SequenceRecord *self, PyObject *args, PyObject *kwargs) {
    PyObject *name = NULL;
    PyObject *sequence = NULL;
    PyObject *qualities = NULL;
    static char * _keywords[] = {"name", "sequence", "qualities", NULL};
    static char * _format = "O!O!|O:SequenceRecord";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, _format, _keywords,
        (PyObject *)&PyUnicode_Type, &name,
        (PyObject *)&PyUnicode_Type, &sequence,
        &qualities))
        return -1;
    if (qualities == Py_None) {
        qualities = NULL;
    }
    if (qualities != NULL) {
        
        if (!PyUnicode_CheckExact(qualities)) {
            PyErr_Format(PyExc_TypeError, "qualities must be of type str, got: %s", Py_TYPE(qualities)->tp_name);
        }
            // Type already checked, can use unsafe macros here.
        if(PyUnicode_GET_LENGTH(sequence) != PyUnicode_GET_LENGTH(qualities)) {
            PyErr_Format(PyExc_ValueError,
                "Size of sequence and qualities do not match: %ld != %ld",
                PyUnicode_GET_LENGTH(sequence), PyUnicode_GET_LENGTH(qualities));
            return -1;
        }
        Py_INCREF(qualities);
    }
    Py_INCREF(name);
    Py_INCREF(sequence);

    self->name = name;
    self->sequence=sequence;
    self->qualities=qualities;
    return 0;
};

PyDoc_STRVAR(BytesSequenceRecord__init____doc__,
"BytesSequenceRecord(name, sequence, qualities)\n"
"--\n"
"\n"
"A sequencing read with read name/id and qualities as bytes objecs\n"
"\n"
"This object only supports FASTQ records." 
"Qualities is a bytes object and it contains the qualities\n"
"encoded as ASCII(qual+33)\n."
"\n"
"  name\n"
"    The read description\n"
"  sequence\n"
"  qualities\n"
"\n"
);
static int
BytesSequenceRecord__init__(SequenceRecord *self, PyObject *args, PyObject *kwargs) {
    PyObject *name = NULL;
    PyObject *sequence = NULL;
    PyObject *qualities = NULL;
    static char * _keywords[] = {"name", "sequence", "qualities", NULL};
    static char * _format = "O!O!O!|:BytesSequenceRecord";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, _format, _keywords, 
        (PyObject *)&PyBytes_Type, &name, 
        (PyObject *)&PyBytes_Type, &sequence, 
        (PyObject *)&PyBytes_Type, &qualities))
        return -1;
    // Type already checked, can use unsafe macros here.
    if(PyBytes_GET_SIZE(sequence) != PyBytes_GET_SIZE(qualities)) {
        PyErr_Format(PyExc_ValueError, 
            "Size of sequence and qualities do not match: %ld != %ld",
            PyBytes_GET_SIZE(sequence), PyBytes_GET_SIZE(qualities));
        return -1;
    }
    Py_INCREF(name);
    Py_INCREF(sequence);
    Py_INCREF(qualities);
    self->name = name; 
    self->sequence=sequence;
    self->qualities=qualities;
    return 0;
};

static PyObject * 
SequenceRecord__repr__(SequenceRecord * self){
    return PyUnicode_FromFormat("SequenceRecord(%R, %R, %R)", 
        self->name, self->sequence, self->qualities);
}

static Py_hash_t 
calculate_sequence_record_hash(SequenceRecord *self){
    Py_hash_t name_hash = PyObject_Hash(self->name);
    Py_hash_t sequence_hash = PyObject_Hash(self->sequence);
    Py_hash_t qualities_hash;
    if (self->qualities != NULL) {
        qualities_hash = PyObject_Hash(self->qualities);
    } else qualities_hash = 0;
    return name_hash ^ sequence_hash ^ qualities_hash;
}

static PyObject *
SequenceRecord__richcompare__(PyObject *self, PyObject *other, int op){
    // This function is extremely generic to allow subtyping, reuse etc.
    if(Py_TYPE(self) != Py_TYPE(other)) {
        PyErr_Format(PyExc_TypeError, 
            "Can only compare objects of %R to objects of the same type. Got: %R.",
            Py_TYPE(self), Py_TYPE(other));
        return NULL;
    }
    if (op == 2){
        return PyBool_FromLong(
            calculate_sequence_record_hash((SequenceRecord *)self) == 
            calculate_sequence_record_hash((SequenceRecord *) other));
    }
    else if (op == 3){
        return PyBool_FromLong(
            calculate_sequence_record_hash((SequenceRecord *)self) != 
            calculate_sequence_record_hash((SequenceRecord *) other));
    }
    else {
        PyErr_Format(PyExc_NotImplementedError, "Only equals and not equals are implemented for %R.", Py_TYPE(self));
        return NULL;
    }
}

static inline PyObject * 
bytes_sequence_to_fastq_record_impl(SequenceRecord *self, int two_headers){
    Py_ssize_t name_length = PyBytes_Size(self->name);
    Py_ssize_t sequence_length = PyBytes_Size(self->sequence);
    Py_ssize_t qualities_length = PyBytes_Size(self->qualities);
    if (name_length == -1 || sequence_length == -1 || qualities_length == -1)
            // Not of type bytes.
            return NULL;
    // Unsafe macros as type is already checked.
    char * name = PyBytes_AS_STRING(self->name);
    char * sequence = PyBytes_AS_STRING(self->sequence);
    char * qualities = PyBytes_AS_STRING(self->qualities);
    return create_fastq_record(name, sequence, qualities,
                               name_length, sequence_length, qualities_length,
                               two_headers);
}

static inline PyObject * 
sequence_to_fastq_record_impl(SequenceRecord *self, int two_headers){
    if (self->qualities == NULL) {
        PyErr_SetString(PyExc_ValueError, 
        "Cannot create FASTQ bytes from a sequence without qualities.");
    }
    Py_ssize_t name_length = PyUnicode_GetLength(self->name);
    Py_ssize_t sequence_length = PyUnicode_GetLength(self->sequence);
    Py_ssize_t qualities_length = PyUnicode_GetLength(self->qualities);
    if (name_length == -1 || sequence_length == -1 || qualities_length == -1)
            // Not of type PyUnicode
            return NULL;
    if (!(
        (PyUnicode_KIND(self->name) == PyUnicode_1BYTE_KIND) && 
        (PyUnicode_KIND(self->sequence)== PyUnicode_1BYTE_KIND) &&
        (PyUnicode_KIND(self->qualities) == PyUnicode_1BYTE_KIND))) {
        PyErr_SetString(PyExc_ValueError, "Name, sequence and qualities must all be valid ASCII strings.");
    }
    // Unsafe macros as type is already checked.
    char * name = (char *)PyUnicode_1BYTE_DATA(self->name);
    char * sequence = (char *)PyUnicode_1BYTE_DATA(self->sequence);
    char * qualities = (char *)PyUnicode_1BYTE_DATA(self->qualities);
    return create_fastq_record(name, sequence, qualities,
                               name_length, sequence_length, qualities_length,
                               two_headers);
}

PyDoc_STRVAR(SequenceRecord_fastq_bytes__doc__,
"Return the entire FASTQ record as bytes which can be written\n"
"into a file.");

#define SEQUENCE_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))SequenceRecord_fastq_bytes, METH_NOARGS, SequenceRecord_fastq_bytes__doc__}

static PyObject *
SequenceRecord_fastq_bytes(SequenceRecord *self, PyObject *NoArgs){
    return sequence_to_fastq_record_impl(self, 0);
}

#define BYTES_SEQUENCE_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))BytesSequenceRecord_fastq_bytes, METH_NOARGS, SequenceRecord_fastq_bytes__doc__}

static PyObject *
BytesSequenceRecord_fastq_bytes(SequenceRecord *self, PyObject *NoArgs){
    return bytes_sequence_to_fastq_record_impl(self, 0);
}
PyDoc_STRVAR(SequenceRecord_fastq_bytes_two_headers__doc__,
"Return this record in FASTQ format as a bytes object where the header\n"
"(after the @) is repeated on the third line.");

#define SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF    \
    {"fastq_bytes_two_headers", (PyCFunction)(void(*)(void))SequenceRecord_fastq_bytes_two_headers, METH_NOARGS, SequenceRecord_fastq_bytes_two_headers__doc__}

static PyObject *
SequenceRecord_fastq_bytes_two_headers(SequenceRecord *self, PyObject *NoArgs){
    return sequence_to_fastq_record_impl(self, 1);
}

#define BYTES_SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF    \
    {"fastq_bytes_two_headers", (PyCFunction)(void(*)(void))BytesSequenceRecord_fastq_bytes_two_headers, METH_NOARGS, SequenceRecord_fastq_bytes_two_headers__doc__}

static PyObject *
BytesSequenceRecord_fastq_bytes_two_headers(SequenceRecord *self, PyObject *NoArgs){
    return bytes_sequence_to_fastq_record_impl(self, 1);
}

static PyMemberDef SequenceRecord_members[] = {
    {"name", T_OBJECT_EX, offsetof(SequenceRecord, name), 0},
    {"sequence", T_OBJECT_EX, offsetof(SequenceRecord, sequence), 0},
    {"qualities", T_OBJECT, offsetof(SequenceRecord, qualities), 0},  // May be None
    {NULL}
};

static PyMemberDef BytesSequenceRecord_members[] = {
    {"name", T_OBJECT_EX, offsetof(SequenceRecord, name), 0},
    {"sequence", T_OBJECT_EX, offsetof(SequenceRecord, sequence), 0},
    {"qualities", T_OBJECT_EX, offsetof(SequenceRecord, qualities), 0},
    {NULL}
};


static PyMethodDef SequenceRecord_methods[] = {
    SEQUENCE_FASTQ_BYTES_METHODDEF,
    SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF,
    {NULL}
};

static PyMethodDef BytesSequenceRecord_methods[] = {
    BYTES_SEQUENCE_FASTQ_BYTES_METHODDEF,
    BYTES_SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF,
    {NULL}
};

static PyTypeObject SequenceRecord_type = {
    .tp_name = "_sequence.SequenceRecord",
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_basicsize = sizeof(SequenceRecord),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)SequenceRecord_dealloc,
    .tp_doc = SequenceRecord__init____doc__,
    .tp_init = (initproc)SequenceRecord__init__,
    .tp_new = PyType_GenericNew,
    .tp_members = SequenceRecord_members,
    .tp_methods = SequenceRecord_methods,
    .tp_repr = (reprfunc)SequenceRecord__repr__,
    .tp_richcompare = SequenceRecord__richcompare__,
};

static PyTypeObject BytesSequenceRecord_type = {
    .tp_name = "_sequence.BytesSequenceRecord",
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_basicsize = sizeof(SequenceRecord),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)SequenceRecord_dealloc,
    .tp_doc = BytesSequenceRecord__init____doc__,
    .tp_init = (initproc)BytesSequenceRecord__init__,
    .tp_new = PyType_GenericNew,
    .tp_members = BytesSequenceRecord_members,
    .tp_methods = BytesSequenceRecord_methods,
    .tp_repr = (reprfunc)SequenceRecord__repr__,
    .tp_richcompare = SequenceRecord__richcompare__,
};


static struct PyModuleDef _sequence_module = {
    PyModuleDef_HEAD_INIT,
    "_sequence",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};


PyMODINIT_FUNC
PyInit__sequence(void)
{
    PyObject *m;

    m = PyModule_Create(&_sequence_module);
    if (m == NULL)
        return NULL;
    PyTypeObject * SequenceRecordType = &SequenceRecord_type;
    PyTypeObject * BytesSequenceRecordType = &BytesSequenceRecord_type;
    if (PyType_Ready(SequenceRecordType) != 0) { 
        return NULL;
    }
    if (PyType_Ready(BytesSequenceRecordType) != 0){
        return NULL;
    }
    Py_INCREF((PyObject *)SequenceRecordType);
    if (PyModule_AddObject(
            m, "SequenceRecord", (PyObject *)SequenceRecordType) != 0) {
        return NULL;
    }
    Py_INCREF((PyObject *)BytesSequenceRecordType);
    if (PyModule_AddObject(
            m, "BytesSequenceRecord", (PyObject *)BytesSequenceRecordType) != 0) {
        return NULL;
    }
    // Add aliases for backwards compatibility
    Py_INCREF((PyObject *)SequenceRecordType);
    if (PyModule_AddObject(
            m, "Sequence", (PyObject *)SequenceRecordType) != 0) {
        return NULL;
    }
    Py_INCREF((PyObject *)BytesSequenceRecordType);
    if (PyModule_AddObject(
            m, "BytesSequence", (PyObject *)BytesSequenceRecordType) != 0) {
        return NULL;
    }
    return m;
}
