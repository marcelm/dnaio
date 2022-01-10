#define PY_SSIZE_T_CLEAN
#include <Python.h>

// This file contains code for creating a Sequence object fast. This bypasses
// the python __init__ and __new__ constructors.
// Given the simplicity of __cinit__ it is simple enough to perform the same
// initialization in using C code.

// Cython generates exactly the same struct under a different name. Unfortunately
// that one cannot be accessed from the cython code.
struct SequenceStruct {
    PyObject_HEAD
    PyObject *name;
    PyObject *sequence;
    PyObject *qualities;
};

static PyObject * new_sequence_object(PyTypeObject *SequenceClass, 
                                      PyObject *name,
                                      PyObject *sequence,
                                      PyObject *qualities)
{
    // PyObject_New ensures the reference count is 1.
    // see https://docs.python.org/3/c-api/allocation.html#c.PyObject_New
    // What happens here is that all the fields of PyObject_HEAD are initialized
    // and that the type pointer is set to point at SequenceClass. Any remaining
    // struct fields not in PyObject_HEAD are not set.
    struct SequenceStruct *new_obj = PyObject_New(struct SequenceStruct, SequenceClass);
    if (new_obj == NULL){
        PyErr_NoMemory();
        return NULL;
    }
    // Increase reference counts so objects do not get deleted while being part of the struct
    // Decreasing the reference counts again will be properly handled in the destructor.
    Py_INCREF(name);
    Py_INCREF(sequence);
    Py_INCREF(qualities);
    new_obj->name = name;
    new_obj->sequence = sequence;
    new_obj->qualities = qualities;
    return (PyObject *)new_obj;
}
