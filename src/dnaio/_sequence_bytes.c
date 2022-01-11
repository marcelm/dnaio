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
SequenceBytes__init__(SequenceBytes *self, PyObject *args, PyObject *kwargs) {
    PyObject *name;
    PyObject *sequence;
    PyObject *qualities;

};

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
    {0, 0},
};

static PyType_Spec Comptype_spec = {
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

    return m;
}
