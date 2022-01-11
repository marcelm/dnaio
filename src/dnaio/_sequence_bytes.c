#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef

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

static PyMemberDef SequenceBytes_members[] = {
    {"name", T_OBJECT_EX, offsetof(SequenceBytes, name), 0},
    {"sequence", T_OBJECT_EX, offsetof(SequenceBytes, sequence), 0},
    {"qualities", T_OBJECT_EX, offsetof(SequenceBytes, qualities), 0},
    {NULL}
};

static PyMethodDef SequenceBytes_methods[] = {
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
    "_sequence_bytes.SequenceBytes",
    sizeof(SequenceBytes_methods),
    0,
    Py_TPFLAGS_DEFAULT,
    SequenceBytes_slots
};


static struct PyModuleDef _sequence_bytes_module = {
    PyModuleDef_HEAD_INIT,
    "_sequence_bytes",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};


PyMODINIT_FUNC
PyInit__sequence_bytes(void)
{
    PyObject *m;

    m = PyModule_Create(&_sequence_bytes_module);
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
