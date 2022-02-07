#define PY_SSIZE_T_CLEAN
#include <Python.h>

typedef struct {
    PyObject_HEAD
    PyObject * name;
    PyObject * sequence;
    PyObject * qualities;
} SequenceRecord;

// NOTE: This function steals references. This is also how CPython tuples work.
// This means that when we create new name, sequence and qualities objects which
// will have their refcount automatically set to 1, we don't need to worry
// about incref or decref. There will be only one reference, and that will be
// owned by the sequence record object.
static PyObject *
new_sequence_record(PyTypeObject *SequenceClass, PyObject *name, PyObject *sequence, PyObject *qualities){
    SequenceRecord *new_obj = PyObject_New(SequenceRecord, SequenceClass);
    if (new_obj == NULL)
        return PyErr_NoMemory();
    new_obj->name = name;
    new_obj->sequence = sequence;
    new_obj->qualities = qualities;
    return (PyObject *)new_obj;
}

static inline PyObject *
create_fastq_record(char * name, char * sequence, char * qualities,
                    Py_ssize_t name_length,
                    Py_ssize_t sequence_length,
                    Py_ssize_t qualities_length,
                    int two_headers) {
    // Total size is name + sequence + qualities + 4 newlines + '+' and an
    // '@' to be put in front of the name.
    Py_ssize_t total_size = name_length + sequence_length + qualities_length + 6;

    if (two_headers)
        // We need space for the name after the +.
        total_size += name_length;

    // This is the canonical way to create an uninitialized bytestring of given size
    PyObject * retval = PyBytes_FromStringAndSize(NULL, total_size);
    if (retval == NULL)
        return PyErr_NoMemory();

    char * retval_ptr = PyBytes_AS_STRING(retval);

    // Write the sequences into the bytestring at the correct positions.
    size_t cursor;
    retval_ptr[0] = '@';
    memcpy(retval_ptr + 1, name, name_length);
    cursor = name_length + 1;
    retval_ptr[cursor] = '\n'; cursor += 1;
    memcpy(retval_ptr + cursor, sequence, sequence_length);
    cursor += sequence_length;
    retval_ptr[cursor] = '\n'; cursor += 1;
    retval_ptr[cursor] = '+'; cursor += 1;
    if (two_headers){
        memcpy(retval_ptr + cursor, name, name_length);
        cursor += name_length;
    }
    retval_ptr[cursor] = '\n'; cursor += 1;
    memcpy(retval_ptr + cursor, qualities, qualities_length);
    cursor += qualities_length;
    retval_ptr[cursor] = '\n';
    return retval;
}