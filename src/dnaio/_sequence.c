#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef
#include "_sequence.h"

typedef struct {
    PyObject_HEAD
    PyObject * name;
    PyObject * sequence;
    PyObject * qualities;
} SequenceBytes;

static void 
SequenceBytes_dealloc(SequenceBytes *self) {
    Py_CLEAR(self->name);
    Py_CLEAR(self->sequence);
    Py_CLEAR(self->qualities);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
new_sequence_bytes(PyTypeObject *SequenceClass, PyObject *name, PyObject *sequence, PyObject *qualities){
    SequenceBytes *new_obj = PyObject_New(SequenceBytes, SequenceClass); 
    if (new_obj == NULL)
        return PyErr_NoMemory();
    Py_INCREF(name);
    Py_INCREF(sequence);
    Py_INCREF(qualities);
    new_obj->name = name;
    new_obj->sequence = sequence;
    new_obj->qualities = qualities;
    return (PyObject *)new_obj;
}

static int
SequenceBytes__init__(SequenceBytes *self, PyObject *args, PyObject *kwargs) {
    PyObject *name = NULL;
    PyObject *sequence = NULL;
    PyObject *qualities = NULL;
    static char * _keywords[] = {"name", "sequence", "qualities", NULL};
    static char * _format = "O!O!O!|:SequenceBytes";
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
SequenceBytes__repr__(SequenceBytes * self){
    return PyUnicode_FromFormat("SequenceBytes(%R, %R, %R)", 
        self->name, self->sequence, self->qualities);
}

static PyObject * 
sequence_bytes_to_fastq_record_impl(SequenceBytes *self, int two_headers){
    Py_ssize_t name_length = PyBytes_Size(self->name);
    Py_ssize_t sequence_length = PyBytes_Size(self->name);
    Py_ssize_t qualities_length = PyBytes_Size(self->name);
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

PyDoc_STRVAR(SequenceBytes_fastq_bytes__doc__,
"Return the entire FASTQ record as bytes which can be written\n"
"into a file.");

#define SEQUENCE_BYTES_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))SequenceBytes_fastq_bytes, METH_NOARGS, SequenceBytes_fastq_bytes__doc__}

static PyObject *
SequenceBytes_fastq_bytes(SequenceBytes *self, PyObject *NoArgs){
    return sequence_bytes_to_fastq_record_impl(self, 0);
}

PyDoc_STRVAR(SequenceBytes_fastq_bytes_two_headers__doc__,
"Return this record in FASTQ format as a bytes object where the header\n"
"(after the @) is repeated on the third line.");

#define SEQUENCE_BYTES_FASTQ_BYTES_TWO_HEADERS_METHODDEF    \
    {"fastq_bytes_two_headers", (PyCFunction)(void(*)(void))SequenceBytes_fastq_bytes_two_headers, METH_NOARGS, SequenceBytes_fastq_bytes_two_headers__doc__}

static PyObject *
SequenceBytes_fastq_bytes_two_headers(SequenceBytes *self, PyObject *NoArgs){
    return sequence_bytes_to_fastq_record_impl(self, 1);
}

static PyMemberDef SequenceBytes_members[] = {
    {"name", T_OBJECT_EX, offsetof(SequenceBytes, name), 0},
    {"sequence", T_OBJECT_EX, offsetof(SequenceBytes, sequence), 0},
    {"qualities", T_OBJECT_EX, offsetof(SequenceBytes, qualities), 0},
    {NULL}
};

static PyMethodDef SequenceBytes_methods[] = {
    SEQUENCE_BYTES_FASTQ_BYTES_METHODDEF,
    SEQUENCE_BYTES_FASTQ_BYTES_TWO_HEADERS_METHODDEF,
    {NULL}
};

static PyType_Slot SequenceBytes_slots[] = {
    {Py_tp_dealloc, SequenceBytes_dealloc},
    {Py_tp_methods, SequenceBytes_methods},
    {Py_tp_members, SequenceBytes_members},
    {Py_tp_doc, "An object containing name, sequence and qualities as bytes objects."},
    {Py_tp_init, SequenceBytes__init__},
    {Py_tp_new, PyType_GenericNew},
    {Py_tp_repr, SequenceBytes__repr__},
    {0, 0},
};

static PyType_Spec SequenceBytes_type_spec = {
    "_sequence.SequenceBytes",
    sizeof(SequenceBytes),
    0,
    Py_TPFLAGS_DEFAULT,
    SequenceBytes_slots
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

    PyTypeObject *SequenceBytesType = (PyTypeObject *)PyType_FromSpec(&SequenceBytes_type_spec);
    if (SequenceBytesType == NULL)
        return NULL;
    
    Py_INCREF(SequenceBytesType);
    if (PyModule_AddObject(m, "SequenceBytes",  (PyObject *)SequenceBytesType) < 0) {
        return NULL;
    }
    return m;
}
