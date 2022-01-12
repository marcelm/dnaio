#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *
new_sequence_bytes(PyTypeObject *SequenceClass, PyObject *name, PyObject *sequence, PyObject *qualities);

#define LINEFEED 10  // \n
#define AT_SYMBOL 64  // @
#define PLUS_SYMBOL 43 // +

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
    retval_ptr[0] = AT_SYMBOL;
    memcpy(retval_ptr + 1, name, name_length);
    cursor = name_length + 1;
    retval_ptr[cursor] = LINEFEED; cursor += 1;
    memcpy(retval_ptr + cursor, sequence, sequence_length);
    cursor += sequence_length;
    retval_ptr[cursor] = LINEFEED; cursor += 1;
    retval_ptr[cursor] = PLUS_SYMBOL; cursor += 1;
    if (two_headers){
        memcpy(retval_ptr + cursor, name, name_length);
        cursor += name_length;
    }
    retval_ptr[cursor] = LINEFEED; cursor += 1;
    memcpy(retval_ptr + cursor, qualities, qualities_length);
    cursor += qualities_length;
    retval_ptr[cursor] = LINEFEED;
    return retval;
}