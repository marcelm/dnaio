# cython: language_level=3, emit_code_comments=False

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING, PyBytes_Check, PyBytes_GET_SIZE
from cpython.unicode cimport PyUnicode_DecodeLatin1, PyUnicode_Check, PyUnicode_GET_LENGTH
from cpython.ref cimport PyObject, Py_INCREF
from libc.string cimport strncmp, memcmp, memcpy, memchr, strcspn
cimport cython

from ._sequence import Sequence, BytesSequence

cdef extern from "_sequence.h":
    object new_sequence_record(type SequenceClass, object name, object sequence, object qualities)

cdef extern from "Python.h":
    unsigned char * PyUnicode_1BYTE_DATA(object o)
    int PyUnicode_KIND(object o)
    int PyUnicode_1BYTE_KIND
    object PyUnicode_New(Py_ssize_t size, Py_UCS4 maxchar)

from typing import Union

from .exceptions import FastqFormatError
from ._util import shorten

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


def fastq_iter(file, sequence_class, Py_ssize_t buffer_size):
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
        bytearray buf = bytearray(buffer_size)
        char[:] buf_view = buf
        char* c_buf = buf
        object name
        object sequence
        object qualities
        Py_ssize_t last_read_position = 0
        Py_ssize_t record_start = 0
        Py_ssize_t bufstart, bufend, name_start, name_end, name_length
        Py_ssize_t sequence_start, sequence_end, sequence_length
        Py_ssize_t second_header_start, second_header_end, second_header_length
        Py_ssize_t qualities_start, qualities_end, qualities_length
        char *name_end_ptr
        char *sequence_end_ptr
        char *second_header_end_ptr
        char *qualities_end_ptr
        bint save_as_bytes = sequence_class is BytesSequence
        bint custom_class = (sequence_class is not Sequence and
                             sequence_class is not BytesSequence)
        Py_ssize_t n_records = 0
        bint extra_newline = False

    if buffer_size < 1:
        raise ValueError("Starting buffer size too small")

    # buf is a byte buffer that is re-used in each iteration. Its layout is:
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

    readinto = file.readinto
    bufstart = 0

    # The input file is processed in chunks that each fit into buf
    while True:
        assert bufstart < len(buf_view)
        bufend = readinto(buf_view[bufstart:]) + bufstart
        if bufstart == bufend:
            # End of file
            if bufstart > 0 and buf_view[bufstart-1] != b'\n':
                # There is still data in the buffer and its last character is
                # not a newline: This is a file that is missing the final
                # newline. Append a newline and continue.
                buf_view[bufstart] = b'\n'
                bufstart += 1
                bufend += 1
                extra_newline = True
            elif last_read_position > record_start:  # Incomplete FASTQ records are present.
                if extra_newline:
                    # Do not report the linefeed that was added by dnaio but
                    # was not present in the original input.
                    last_read_position -= 1
                lines = buf[record_start:last_read_position].count(b'\n')
                raise FastqFormatError(
                    'Premature end of file encountered. The incomplete final record was: '
                    '{!r}'.format(
                        shorten(buf[record_start:last_read_position].decode('latin-1'),
                                500)),
                    line=n_records * 4 + lines)
            else:  # EOF Reached. Stop iterating.
                return

        # Parse all complete FASTQ records in this chunk
        record_start = 0
        while True:
            ### Check for a complete record (i.e 4 newlines are present)
            # Use libc memchr, this optimizes looking for characters by
            # using 64-bit integers. See:
            # https://sourceware.org/git/?p=glibc.git;a=blob_plain;f=string/memchr.c;hb=HEAD
            # void *memchr(const void *str, int c, size_t n)
            name_end_ptr = <char *>memchr(c_buf + record_start, b'\n', <size_t>(bufend - record_start))
            if name_end_ptr == NULL:
                break
            # bufend - sequence_start is always nonnegative:
            # - name_end is at most bufend - 1
            # - thus sequence_start is at most bufend
            name_end = name_end_ptr - c_buf
            sequence_start = name_end + 1
            sequence_end_ptr = <char *>memchr(c_buf + sequence_start, b'\n', <size_t>(bufend - sequence_start))
            if sequence_end_ptr == NULL:
                break
            sequence_end = sequence_end_ptr - c_buf
            second_header_start = sequence_end + 1
            second_header_end_ptr = <char *>memchr(c_buf + second_header_start, b'\n', <size_t>(bufend - second_header_start))
            if second_header_end_ptr == NULL:
                break
            second_header_end = second_header_end_ptr - c_buf
            qualities_start = second_header_end + 1
            qualities_end_ptr = <char *>memchr(c_buf + qualities_start, b'\n', <size_t>(bufend - qualities_start))
            if qualities_end_ptr == NULL:
                break
            qualities_end = qualities_end_ptr - c_buf

            if c_buf[record_start] != b'@':
                raise FastqFormatError("Line expected to "
                    "start with '@', but found {!r}".format(chr(c_buf[record_start])),
                    line=n_records * 4)
            if c_buf[second_header_start] != b'+':
                raise FastqFormatError("Line expected to "
                    "start with '+', but found {!r}".format(chr(c_buf[second_header_start])),
                    line=n_records * 4 + 2)

            name_start = record_start + 1  # Skip @
            second_header_start += 1  # Skip +
            name_length = name_end - name_start
            sequence_length = sequence_end - sequence_start
            second_header_length = second_header_end - second_header_start
            qualities_length = qualities_end - qualities_start

            # Check for \r\n line-endings and compensate
            if c_buf[name_end - 1] == b'\r':
                name_length -= 1
            if c_buf[sequence_end - 1] == b'\r':
                sequence_length -= 1
            if c_buf[second_header_end - 1] == b'\r':
                second_header_length -= 1
            if c_buf[qualities_end - 1] == b'\r':
                qualities_length -= 1

            if second_header_length:  # should be 0 when only + is present
                if (name_length != second_header_length or
                        strncmp(c_buf+second_header_start,
                            c_buf + name_start, second_header_length) != 0):
                    raise FastqFormatError(
                        "Sequence descriptions don't match ('{}' != '{}').\n"
                        "The second sequence description must be either "
                        "empty or equal to the first description.".format(
                            c_buf[name_start:name_end].decode('latin-1'),
                            c_buf[second_header_start:second_header_end]
                            .decode('latin-1')), line=n_records * 4 + 2)

            if qualities_length != sequence_length:
                raise FastqFormatError(
                    "Length of sequence and qualities differ", line=n_records * 4 + 3)

            if n_records == 0:
                yield bool(second_header_length)  # first yielded value is special

            if save_as_bytes:
                name = PyBytes_FromStringAndSize(c_buf + name_start, name_length)
                sequence = PyBytes_FromStringAndSize(c_buf + sequence_start, sequence_length)
                qualities = PyBytes_FromStringAndSize(c_buf + qualities_start, qualities_length)
            else:
                # Constructing objects with PyUnicode_New and memcpy bypasses some of
                # the checks otherwise done when using PyUnicode_DecodeLatin1 or similar
                name = PyUnicode_New(name_length, 255)
                sequence = PyUnicode_New(sequence_length, 255)
                qualities = PyUnicode_New(qualities_length, 255)
                if <PyObject*>name == NULL or <PyObject*>sequence == NULL or <PyObject*>qualities == NULL:
                    raise MemoryError()
                memcpy(PyUnicode_1BYTE_DATA(name), c_buf + name_start, name_length)
                memcpy(PyUnicode_1BYTE_DATA(sequence), c_buf + sequence_start, sequence_length)
                memcpy(PyUnicode_1BYTE_DATA(qualities), c_buf + qualities_start, qualities_length)
                
            if custom_class:
                yield sequence_class(name, sequence, qualities)
            else:
                Py_INCREF(name); Py_INCREF(sequence); Py_INCREF(qualities)
                yield new_sequence_record(sequence_class, name, sequence, qualities)

            ### Advance record to next position
            n_records += 1
            record_start = qualities_end + 1
        # bufend reached
        last_read_position = bufend
        if record_start == 0 and bufend == len(buf):
            # buffer too small, double it
            buffer_size *= 2
            prev_buf = buf
            buf = bytearray(buffer_size)
            buf[0:bufend] = prev_buf
            del prev_buf
            bufstart = bufend
            buf_view = buf
            c_buf = buf
        else:
            bufstart = bufend - record_start
            buf[0:bufstart] = buf[record_start:bufend]


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
