# cython: language_level=3, emit_code_comments=False

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING
from libc.string cimport strncmp, memcmp, memcpy, memchr, strcspn
from cpython.unicode cimport PyUnicode_GET_LENGTH

from .exceptions import FastqFormatError
from ._util import shorten


ctypedef public class Sequence[type SequenceType, object SequenceStruct]:
    """
    A sequencing read with read name/id and (optional) qualities

    If qualities are available, they are as
    For a Sequence a FASTA file
    record containing a read in a FASTA or FASTQ file. For FASTA, the qualities attribute
    is None. For FASTQ, qualities is a string and it contains the qualities
    encoded as ASCII(qual+33).

    Attributes:
      name (str): The read description
      sequence (str):
      qualities (str):
    """
    cdef:
        public str name
        public str sequence
        public str qualities

    def __cinit__(self, str name, str sequence, str qualities=None):
        """Set qualities to None if there are no quality values"""
        self.name = name
        self.sequence = sequence
        self.qualities = qualities

    def __init__(self, str name, str sequence, str qualities = None):
        # __cinit__ is called first and sets all the variables.
        if qualities is not None and len(qualities) != len(sequence):
            rname = shorten(name)
            raise ValueError("In read named {!r}: length of quality sequence "
                             "({}) and length of read ({}) do not match".format(
                rname, len(qualities), len(sequence)))

    def __getitem__(self, key):
        """
        Slice this Sequence. If the qualities attribute is not None, it is
        sliced accordingly. The read name is copied unchanged.

        Returns:
          A new Sequence object with a sliced sequence.
        """
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.qualities is not None else None)

    def __repr__(self):
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={!r}'.format(shorten(self.qualities))
        return '<Sequence(name={!r}, sequence={!r}{})>'.format(
            shorten(self.name), shorten(self.sequence), qstr)

    def __len__(self):
        """
        Returns:
           The number of characters in this sequence
        """
        return len(self.sequence)

    def __richcmp__(self, other, int op):
        if 2 <= op <= 3:
            eq = self.name == other.name and \
                self.sequence == other.sequence and \
                self.qualities == other.qualities
            if op == 2:
                return eq
            else:
                return not eq
        else:
            raise NotImplementedError()

    def __reduce__(self):
        return (Sequence, (self.name, self.sequence, self.qualities))

    def qualities_as_bytes(self):
        """Return the qualities as a bytes object.

        This is a faster version of qualities.encode('ascii')."""
        return self.qualities.encode('ascii')

    def fastq_bytes(self):
        """Return the entire FASTQ record as bytes which can be written
        into a file."""
        # Convert to ASCII bytes sequences first as these have a one-to-one
        # relation between size and number of bytes
        # Unlike decoding, ascii is not slower than latin-1. This is because
        # CPython performs a call to PyUnicodeCheck on both occassions. This
        # determines the type of the Unicode object. In fact, the ascii encode
        # is slightly faster because the check for PyASCIIObject is performed
        # first.
        cdef bytes name = self.name.encode('ascii')
        cdef bytes sequence = self.sequence.encode('ascii')
        cdef bytes qualities = self.qualities.encode('ascii')
        cdef Py_ssize_t name_length = len(name)
        cdef Py_ssize_t sequence_length = len(sequence)
        cdef Py_ssize_t qualities_length = len(qualities)

        # Since Cython will generate code above that is a 100% sure to generate
        # bytes objects, we can call Python C-API functions that don't perform
        # checks on the object.
        cdef char * name_ptr = PyBytes_AS_STRING(name)
        cdef char * sequence_ptr = PyBytes_AS_STRING(sequence)
        cdef char * qualities_ptr = PyBytes_AS_STRING(qualities)

        # Total size is name + sequence + qualities + 4 newlines + '+' and an
        # '@' to be put in front of the name.
        cdef Py_ssize_t total_size = name_length + sequence_length + qualities_length + 6

        # This is the canonical way to create an uninitialized bytestring of given size
        cdef bytes retval = PyBytes_FromStringAndSize(NULL, total_size)
        cdef char * retval_ptr = PyBytes_AS_STRING(retval)

        # Write the sequences into the bytestring at the correct positions.
        cdef Py_ssize_t cursor
        retval_ptr[0] = b"@"
        memcpy(retval_ptr + 1, name_ptr, name_length)
        cursor = name_length + 1
        retval_ptr[cursor] = b"\n"; cursor += 1
        memcpy(retval_ptr + cursor, sequence_ptr, sequence_length)
        cursor += sequence_length
        retval_ptr[cursor] = b"\n"; cursor += 1
        retval_ptr[cursor] = b"+"; cursor += 1
        retval_ptr[cursor] = b"\n"; cursor += 1
        memcpy(retval_ptr + cursor, qualities_ptr, qualities_length)
        cursor += qualities_length
        retval_ptr[cursor] = b"\n"
        return retval

    def fastq_bytes_two_headers(self):
        """
        Return this record in FASTQ format as a bytes object where the header (after the @) is
        repeated on the third line.
        """
        return f"@{self.name}\n{self.sequence}\n+{self.name}\n{self.qualities}\n".encode("ascii")

