# cython: language_level=3, emit_code_comments=False

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING, PyBytes_Check, PyBytes_GET_SIZE
from cpython.unicode cimport PyUnicode_DecodeLatin1, PyUnicode_Check, PyUnicode_GET_LENGTH
from cpython.ref cimport PyObject
from libc.string cimport strncmp, memcmp, memcpy, memchr, strcspn
cimport cython

cdef extern from "Python.h":
    unsigned char * PyUnicode_1BYTE_DATA(object o)
    int PyUnicode_KIND(object o)
    int PyUnicode_1BYTE_KIND
    object PyUnicode_New(Py_ssize_t size, Py_UCS4 maxchar)

from typing import Union

from .exceptions import FastqFormatError
from ._util import shorten


cdef class Sequence:
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

    def fastq_bytes(self, two_headers = False):
        """Return the entire FASTQ record as bytes which can be written
        into a file.

        Optionally the header (after the @) can be repeated on the third line
        (after the +), when two_headers is enabled."""
        cdef:
            char * name
            char * sequence
            char * qualities
            Py_ssize_t name_length
            Py_ssize_t sequence_length
            Py_ssize_t qualities_length

        if PyUnicode_KIND(self.name) == PyUnicode_1BYTE_KIND:
            name = <char *>PyUnicode_1BYTE_DATA(self.name)
            name_length = <size_t>PyUnicode_GET_LENGTH(self.name)
        else:
            # Allow non-ASCII in name
            name_bytes = self.name.encode('latin-1')
            name = PyBytes_AS_STRING(name_bytes)
            name_length = PyBytes_GET_SIZE(name_bytes)
        if PyUnicode_KIND(self.sequence) == PyUnicode_1BYTE_KIND:
            sequence = <char *>PyUnicode_1BYTE_DATA(self.sequence)
            sequence_length = <size_t>PyUnicode_GET_LENGTH(self.sequence)
        else:
            # Don't allow non-ASCII in sequence and qualities
            sequence_bytes = self.sequence.encode('ascii')
            sequence = PyBytes_AS_STRING(sequence_bytes)
            sequence_length = PyBytes_GET_SIZE(sequence_bytes)
        if PyUnicode_KIND(self.qualities) == PyUnicode_1BYTE_KIND:
            qualities = <char *>PyUnicode_1BYTE_DATA(self.qualities)
            qualities_length = <size_t>PyUnicode_GET_LENGTH(self.qualities)
        else:
            qualities_bytes = self.qualities.encode('ascii')
            qualities = PyBytes_AS_STRING(qualities_bytes)
            qualities_length = PyBytes_GET_SIZE(qualities_bytes)
        return create_fastq_record(name, sequence, qualities,
                                   name_length, sequence_length, qualities_length,
                                   two_headers)

    def fastq_bytes_two_headers(self):
        """
        Return this record in FASTQ format as a bytes object where the header (after the @) is
        repeated on the third line.
        """
        return self.fastq_bytes(two_headers=True)


cdef class BytesSequence:
    cdef:
        public bytes name
        public bytes sequence
        public bytes qualities

    def __cinit__(self, bytes name, bytes sequence, bytes qualities):
        """Set qualities to None if there are no quality values"""
        self.name = name
        self.sequence = sequence
        self.qualities = qualities

    def __init__(self, bytes name, bytes sequence, bytes qualities):
        # __cinit__ is called first and sets all the variables.
        if len(qualities) != len(sequence):
            rname = shorten(name)
            raise ValueError("In read named {!r}: length of quality sequence "
                             "({}) and length of read ({}) do not match".format(
                rname, len(qualities), len(sequence)))
    
    def __repr__(self):
        return '<BytesSequence(name={!r}, sequence={!r}, qualities={!r})>'.format(
            shorten(self.name), shorten(self.sequence), shorten(self.qualities))

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

    def fastq_bytes(self, two_headers=False):
        name = PyBytes_AS_STRING(self.name)
        name_length = PyBytes_GET_SIZE(self.name)
        sequence = PyBytes_AS_STRING(self.sequence)
        sequence_length = PyBytes_GET_SIZE(self.sequence)
        qualities = PyBytes_AS_STRING(self.qualities)
        qualities_length = PyBytes_GET_SIZE(self.qualities)
        return create_fastq_record(name, sequence, qualities,
                                   name_length, sequence_length, qualities_length,
                                   two_headers)

    def fastq_bytes_two_headers(self):
        """
        Return this record in FASTQ format as a bytes object where the header (after the @) is
        repeated on the third line.
        """
        return self.fastq_bytes(two_headers=True)


cdef bytes create_fastq_record(char * name, char * sequence, char * qualities,
                               Py_ssize_t name_length,
                               Py_ssize_t sequence_length,
                               Py_ssize_t qualities_length,
                               bint two_headers = False):
        # Total size is name + sequence + qualities + 4 newlines + '+' and an
        # '@' to be put in front of the name.
        cdef Py_ssize_t total_size = name_length + sequence_length + qualities_length + 6

        if two_headers:
            # We need space for the name after the +.
            total_size += name_length

        # This is the canonical way to create an uninitialized bytestring of given size
        cdef bytes retval = PyBytes_FromStringAndSize(NULL, total_size)
        cdef char * retval_ptr = PyBytes_AS_STRING(retval)

        # Write the sequences into the bytestring at the correct positions.
        cdef size_t cursor
        retval_ptr[0] = b"@"
        memcpy(retval_ptr + 1, name, name_length)
        cursor = name_length + 1
        retval_ptr[cursor] = b"\n"; cursor += 1
        memcpy(retval_ptr + cursor, sequence, sequence_length)
        cursor += sequence_length
        retval_ptr[cursor] = b"\n"; cursor += 1
        retval_ptr[cursor] = b"+"; cursor += 1
        if two_headers:
            memcpy(retval_ptr + cursor, name, name_length)
            cursor += name_length
        retval_ptr[cursor] = b"\n"; cursor += 1
        memcpy(retval_ptr + cursor, qualities, qualities_length)
        cursor += qualities_length
        retval_ptr[cursor] = b"\n"
        return retval

# It would be nice to be able to have the first parameter be an
# unsigned char[:] (memory view), but this fails with a BufferError
# when a bytes object is passed in.
# See <https://stackoverflow.com/questions/28203670/>

ctypedef fused bytes_or_bytearray:
    bytes
    bytearray


def paired_fastq_heads(bytes_or_bytearray buf1, bytes_or_bytearray buf2, Py_ssize_t end1, Py_ssize_t end2):
    """
    Skip forward in the two buffers by multiples of four lines.

    Return a tuple (length1, length2) such that buf1[:length1] and
    buf2[:length2] contain the same number of lines (where the
    line number is divisible by four).
    """
    cdef:
        Py_ssize_t pos1 = 0, pos2 = 0
        Py_ssize_t linebreaks = 0
        unsigned char* data1 = buf1
        unsigned char* data2 = buf2
        Py_ssize_t record_start1 = 0
        Py_ssize_t record_start2 = 0

    while True:
        while pos1 < end1 and data1[pos1] != b'\n':
            pos1 += 1
        if pos1 == end1:
            break
        pos1 += 1
        while pos2 < end2 and data2[pos2] != b'\n':
            pos2 += 1
        if pos2 == end2:
            break
        pos2 += 1
        linebreaks += 1
        if linebreaks == 4:
            linebreaks = 0
            record_start1 = pos1
            record_start2 = pos2

    # Hit the end of the data block
    return record_start1, record_start2


cdef class fastq_iter:
    """
    Parse a FASTQ file and yield Sequence objects

    The *first value* that the generator yields is a boolean indicating whether
    the first record in the FASTQ has a repeated header (in the third row
    after the ``+``).

    file -- a file-like object, opened in binary mode (it must have a readinto
    method)

    buffer_size -- size of the initial buffer. This is automatically grown
        if a FASTQ record is encountered that does not fit.
    """
    cdef:
        Py_ssize_t buffer_size
        bytearray buf
        char[:] buf_view
        char *c_buf
        type sequence_class
        bint save_as_bytes
        bint custom_class
        bint extra_newline
        bint yielded_two_headers
        bint eof
        object file
        Py_ssize_t bufend
        Py_ssize_t record_start
    cdef readonly Py_ssize_t n_records

    def __cinit__(self, file, sequence_class, Py_ssize_t buffer_size):
        self.buffer_size = buffer_size
        self.buf = bytearray(buffer_size)
        self.buf_view = self.buf
        self.c_buf = self.buf
        self.sequence_class = sequence_class
        self.save_as_bytes = sequence_class is BytesSequence
        self.custom_class = (sequence_class is not Sequence and
                             sequence_class is not BytesSequence)
        self.n_records = 0
        self.extra_newline = False
        self.yielded_two_headers = False
        self.eof = False
        self.bufend = 0
        self.record_start = 0
        self.file = file
        if buffer_size < 1:
            raise ValueError("Starting buffer size too small")

    cdef _update_buffer(self):
        # self.buf is a byte buffer that is re-used in each iteration. Its layout is:
        #
        # |-- complete records --|
        # +---+------------------+---------+-------+
        # |   |                  |         |       |
        # +---+------------------+---------+-------+
        # ^   ^                  ^         ^       ^
        # 0   bufstart           end       bufend  len(buf)
        #
        # buf[0:bufstart] is the 'leftover' data that could not be processed
        # in the previous iteration because it contained an incomplete
        # FASTQ record.

        cdef Py_ssize_t last_read_position = self.bufend
        cdef Py_ssize_t bufstart
        if self.record_start == 0 and self.bufend == len(self.buf):
            # buffer too small, double it
            self.buffer_size *= 2
            prev_buf = self.buf
            self.buf = bytearray(self.buffer_size)
            self.buf[0:self.bufend] = prev_buf
            del prev_buf
            bufstart = self.bufend
            self.buf_view = self.buf
            self.c_buf = self.buf
        else:
            bufstart = self.bufend - self.record_start
            self.buf[0:bufstart] = self.buf[self.record_start:self.bufend]
        assert bufstart < len(self.buf_view)
        self.bufend = self.file.readinto(self.buf_view[bufstart:]) + bufstart
        if bufstart == self.bufend:
            # End of file
            if bufstart > 0 and self.buf_view[bufstart-1] != b'\n':
                # There is still data in the buffer and its last character is
                # not a newline: This is a file that is missing the final
                # newline. Append a newline and continue.
                self.buf_view[bufstart] = b'\n'
                bufstart += 1
                self.bufend += 1
                self.extra_newline = True
            elif last_read_position > self.record_start:  # Incomplete FASTQ records are present.
                if self.extra_newline:
                    # Do not report the linefeed that was added by dnaio but
                    # was not present in the original input.
                    last_read_position -= 1
                lines = self.buf[self.record_start:last_read_position].count(b'\n')
                raise FastqFormatError(
                    'Premature end of file encountered. The incomplete final record was: '
                    '{!r}'.format(
                        shorten(self.buf[self.record_start:last_read_position].decode('latin-1'),
                                500)),
                    line=self.n_records * 4 + lines)
            else:  # EOF Reached. Stop iterating.
                self.eof = True
        self.record_start = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef:
            object ret_val
            Py_ssize_t bufstart, bufend, name_start, name_end, name_length
            Py_ssize_t sequence_start, sequence_end, sequence_length
            Py_ssize_t second_header_start, second_header_end, second_header_length
            Py_ssize_t qualities_start, qualities_end, qualities_length
            char *name_end_ptr
            char *sequence_end_ptr
            char *second_header_end_ptr
            char *qualities_end_ptr
        while True:
            if self.eof:
                raise StopIteration()
            ### Check for a complete record (i.e 4 newlines are present)
            # Use libc memchr, this optimizes looking for characters by
            # using 64-bit integers. See:
            # https://sourceware.org/git/?p=glibc.git;a=blob_plain;f=string/memchr.c;hb=HEAD
            # void *memchr(const void *str, int c, size_t n)
            name_end_ptr = <char *>memchr(self.c_buf + self.record_start, b'\n', <size_t>(self.bufend - self.record_start))
            if name_end_ptr == NULL:
                self._update_buffer()
                continue
            # bufend - sequence_start is always nonnegative:
            # - name_end is at most bufend - 1
            # - thus sequence_start is at most bufend
            name_end = name_end_ptr - self.c_buf
            sequence_start = name_end + 1
            sequence_end_ptr = <char *>memchr(self.c_buf + sequence_start, b'\n', <size_t>(self.bufend - sequence_start))
            if sequence_end_ptr == NULL:
                self._update_buffer()
                continue
            sequence_end = sequence_end_ptr - self.c_buf
            second_header_start = sequence_end + 1
            second_header_end_ptr = <char *>memchr(self.c_buf + second_header_start, b'\n', <size_t>(self.bufend - second_header_start))
            if second_header_end_ptr == NULL:
                self._update_buffer()
                continue
            second_header_end = second_header_end_ptr - self.c_buf
            qualities_start = second_header_end + 1
            qualities_end_ptr = <char *>memchr(self.c_buf + qualities_start, b'\n', <size_t>(self.bufend - qualities_start))
            if qualities_end_ptr == NULL:
                self._update_buffer()
                continue
            qualities_end = qualities_end_ptr - self.c_buf

            if self.c_buf[self.record_start] != b'@':
                raise FastqFormatError("Line expected to "
                    "start with '@', but found {!r}".format(chr(self.c_buf[self.record_start])),
                    line=self.n_records * 4)
            if self.c_buf[second_header_start] != b'+':
                raise FastqFormatError("Line expected to "
                    "start with '+', but found {!r}".format(chr(self.c_buf[second_header_start])),
                    line=self.n_records * 4 + 2)

            name_start = self.record_start + 1  # Skip @
            second_header_start += 1  # Skip +
            name_length = name_end - name_start
            sequence_length = sequence_end - sequence_start
            second_header_length = second_header_end - second_header_start
            qualities_length = qualities_end - qualities_start

            # Check for \r\n line-endings and compensate
            if self.c_buf[name_end - 1] == b'\r':
                name_length -= 1
            if self.c_buf[sequence_end - 1] == b'\r':
                sequence_length -= 1
            if self.c_buf[second_header_end - 1] == b'\r':
                second_header_length -= 1
            if self.c_buf[qualities_end - 1] == b'\r':
                qualities_length -= 1

            if second_header_length:  # should be 0 when only + is present
                if (name_length != second_header_length or
                        strncmp(self.c_buf+second_header_start,
                            self.c_buf + name_start, second_header_length) != 0):
                    raise FastqFormatError(
                        "Sequence descriptions don't match ('{}' != '{}').\n"
                        "The second sequence description must be either "
                        "empty or equal to the first description.".format(
                            self.c_buf[name_start:name_end].decode('latin-1'),
                            self.c_buf[second_header_start:second_header_end]
                            .decode('latin-1')), line=self.n_records * 4 + 2)

            if qualities_length != sequence_length:
                raise FastqFormatError(
                    "Length of sequence and qualities differ", line=self.n_records * 4 + 3)

            if self.n_records == 0 and not self.yielded_two_headers:
                self.yielded_two_headers = True
                return bool(second_header_length)  # first yielded value is special

            if self.save_as_bytes:
                name = PyBytes_FromStringAndSize(self.c_buf + name_start, name_length)
                sequence = PyBytes_FromStringAndSize(self.c_buf + sequence_start, sequence_length)
                qualities = PyBytes_FromStringAndSize(self.c_buf + qualities_start, qualities_length)
                ret_val = BytesSequence.__new__(BytesSequence, name, sequence, qualities)
            else:
                # Constructing objects with PyUnicode_New and memcpy bypasses some of
                # the checks otherwise done when using PyUnicode_DecodeLatin1 or similar
                name = PyUnicode_New(name_length, 255)
                sequence = PyUnicode_New(sequence_length, 255)
                qualities = PyUnicode_New(qualities_length, 255)
                if <PyObject*>name == NULL or <PyObject*>sequence == NULL or <PyObject*>qualities == NULL:
                    raise MemoryError()
                memcpy(PyUnicode_1BYTE_DATA(name), self.c_buf + name_start, name_length)
                memcpy(PyUnicode_1BYTE_DATA(sequence), self.c_buf + sequence_start, sequence_length)
                memcpy(PyUnicode_1BYTE_DATA(qualities), self.c_buf + qualities_start, qualities_length)
                
                if self.custom_class:
                    ret_val = self.sequence_class(name, sequence, qualities)
                else:
                    ret_val = Sequence.__new__(Sequence, name, sequence, qualities)

            ### Advance record to next position
            self.n_records += 1
            self.record_start = qualities_end + 1
            return ret_val


def record_names_match(header1: str, header2: str):
    """
    Check whether the sequence record ids id1 and id2 are compatible, ignoring a
    suffix of '1', '2' or '3'. This exception allows to check some old
    paired-end reads that have IDs ending in '/1' and '/2'. Also, the
    fastq-dump tool (used for converting SRA files to FASTQ) appends '.1', '.2'
    and sometimes '.3' to paired-end reads if option -I is used.
    """
    cdef:
        char * header1_chars = NULL
        char * header2_chars = NULL
        size_t header1_length
    if PyUnicode_Check(header1):
        if PyUnicode_KIND(header1) == PyUnicode_1BYTE_KIND:
            header1_chars = <char *>PyUnicode_1BYTE_DATA(header1)
            header1_length = <size_t> PyUnicode_GET_LENGTH(header1)
        else:
            header1 = header1.encode('latin1')
            header1_chars = PyBytes_AS_STRING(header1)
            header1_length = PyBytes_GET_SIZE(header1)
    else:
        raise TypeError(f"Header 1 is the wrong type. Expected bytes or string, "
                        f"got: {type(header1)}")

    if PyUnicode_Check(header2):
        if PyUnicode_KIND(header2) == PyUnicode_1BYTE_KIND:
            header2_chars = <char *>PyUnicode_1BYTE_DATA(header2)
        else:
            header2 = header2.encode('latin1')
            header2_chars = PyBytes_AS_STRING(header2)
    else:
        raise TypeError(f"Header 2 is the wrong type. Expected bytes or string, "
                        f"got: {type(header2)}")

    return record_ids_match(header1_chars, header2_chars, header1_length)


def record_names_match_bytes(header1: bytes, header2: bytes):
    if not (PyBytes_Check(header1) and PyBytes_Check(header2)):
        raise TypeError("Header1 and header2 should both be bytes objects. "
                        "Got {} and {}".format(type(header1), type(header2)))
    return record_ids_match(PyBytes_AS_STRING(header1),
                            PyBytes_AS_STRING(header2),
                            PyBytes_GET_SIZE(header1))

cdef bint record_ids_match(char *header1, char *header2, size_t header1_length):
    """
    Check whether the ASCII-encoded IDs match. Only header1_length is needed.
    """
    # Only the read ID is of interest.
    # Find the first tab or space, if not present, strcspn will return the
    # position of the terminating NULL byte. (I.e. the length).
    # Header1 is not searched because we can reuse the end of ID position of
    # header2 as header1's ID should end at the same position.
    cdef size_t id2_length = strcspn(header2, b' \t')

    if header1_length < id2_length:
        return False

    cdef char end = header1[id2_length]
    if end != b'\000' and end != b' ' and end != b'\t':
        return False

    # Check if the IDs end with 1, 2 or 3. This is the read pair number
    # which should not be included in the comparison.
    cdef bint id1endswithnumber = b'1' <= header1[id2_length - 1] <= b'3'
    cdef bint id2endswithnumber = b'1' <= header2[id2_length - 1] <= b'3'
    if id1endswithnumber and id2endswithnumber:
        id2_length -= 1

    # Compare the strings up to the ID end position.
    return memcmp(<void *>header1, <void *>header2, id2_length) == 0
