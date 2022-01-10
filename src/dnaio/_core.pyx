# cython: language_level=3, emit_code_comments=False

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING
from libc.string cimport strncmp, memcmp, memcpy, memchr, strcspn
from cpython.unicode cimport PyUnicode_GET_LENGTH
cimport cython

cdef extern from *:
    unsigned char * PyUnicode_1BYTE_DATA(object o)
    int PyUnicode_KIND(object o)
    int PyUnicode_1BYTE_KIND

from .exceptions import FastqFormatError
from ._util import shorten
from ._sequence import Sequence

cdef extern from "_new_sequence.h":
    object new_sequence_object(type SequenceClass, object name, object sequence, object qualities)


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
        int endskip
        str name
        char* name_encoded
        Py_ssize_t bufstart, bufend, pos, record_start, sequence_start
        Py_ssize_t second_header_start, sequence_length, qualities_start
        Py_ssize_t second_header_length, name_length
        bint custom_class = sequence_class is not Sequence
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
            else:
                break

        # Parse all complete FASTQ records in this chunk
        pos = 0
        record_start = 0
        while True:
            # Parse the name (line 0)
            if c_buf[pos] != b'@':
                raise FastqFormatError("Line expected to "
                    "start with '@', but found {!r}".format(chr(c_buf[pos])),
                    line=n_records * 4)
            pos += 1
            while pos < bufend and c_buf[pos] != b'\n':
                pos += 1
            if pos == bufend:
                break
            endskip = 1 if c_buf[pos-1] == b'\r' else 0
            name_length = pos - endskip - record_start - 1
            name_encoded = c_buf + record_start + 1
            # .decode('latin-1') is 50% faster than .decode('ascii')
            # This is because PyUnicode_DecodeLatin1 is an alias for
            # _PyUnicode_FromUCS1. Which directly copies the bytes into a
            # string object. No operations are taking place. With
            # PyUnicode_DecodeASCII, all characters are checked whether they
            # exceed 128.
            name = c_buf[record_start+1:pos-endskip].decode('latin-1')

            pos += 1

            # Parse the sequence (line 1)
            sequence_start = pos
            while pos < bufend and c_buf[pos] != b'\n':
                pos += 1
            if pos == bufend:
                break
            endskip = 1 if c_buf[pos-1] == b'\r' else 0
            sequence = c_buf[sequence_start:pos-endskip].decode('latin-1')
            sequence_length = pos - endskip - sequence_start
            pos += 1

            # Parse second header (line 2)
            second_header_start = pos
            if pos == bufend:
                break
            if c_buf[pos] != b'+':
                raise FastqFormatError("Line expected to "
                    "start with '+', but found {!r}".format(chr(c_buf[pos])),
                    line=n_records * 4 + 2)
            pos += 1  # skip over the '+'
            while pos < bufend and c_buf[pos] != b'\n':
                pos += 1
            if pos == bufend:
                break
            endskip = 1 if c_buf[pos-1] == b'\r' else 0
            second_header_length = pos - endskip - second_header_start - 1
            if second_header_length == 0:
                second_header = False
            else:
                if (name_length != second_header_length or
                        strncmp(c_buf+second_header_start+1,
                            name_encoded, second_header_length) != 0):
                    raise FastqFormatError(
                        "Sequence descriptions don't match ('{}' != '{}').\n"
                        "The second sequence description must be either "
                        "empty or equal to the first description.".format(
                            name_encoded[:name_length].decode('latin-1'),
                            c_buf[second_header_start+1:pos-endskip]
                            .decode('latin-1')), line=n_records * 4 + 2)
                second_header = True
            pos += 1

            # Parse qualities (line 3)
            qualities_start = pos
            while pos < bufend and c_buf[pos] != b'\n':
                pos += 1
            if pos == bufend:
                break
            endskip = 1 if c_buf[pos-1] == b'\r' else 0
            qualities = c_buf[qualities_start:pos-endskip].decode('latin-1')
            if pos - endskip - qualities_start != sequence_length:
                raise FastqFormatError("Length of sequence and "
                    "qualities differ", line=n_records * 4 + 3)
            pos += 1
            if n_records == 0:
                yield second_header  # first yielded value is special
            if custom_class:
                yield sequence_class(name, sequence, qualities)
            else:
                yield new_sequence_object(Sequence, name, sequence, qualities)
            n_records += 1
            record_start = pos
            if pos == bufend:
                break
        if pos == bufend:
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
    if pos > record_start:
        if extra_newline:
            pos -= 1
        lines = buf[record_start:pos].count(b'\n')
        raise FastqFormatError(
            'Premature end of file encountered. The incomplete final record was: '
            '{!r}'.format(shorten(buf[record_start:pos].decode('latin-1'), 500)),
            line=n_records * 4 + lines)


def record_names_match(header1: str, header2: str):
    """
    Check whether the sequence record ids id1 and id2 are compatible, ignoring a
    suffix of '1', '2' or '3'. This exception allows to check some old
    paired-end reads that have IDs ending in '/1' and '/2'. Also, the
    fastq-dump tool (used for converting SRA files to FASTQ) appends '.1', '.2'
    and sometimes '.3' to paired-end reads if option -I is used.
    """
    if (
        PyUnicode_KIND(header1) != PyUnicode_1BYTE_KIND or
        PyUnicode_KIND(header2) != PyUnicode_1BYTE_KIND
    ):
        # Fall back to slower code path.
        name1 = header1.split(maxsplit=1)[0]
        name2 = header2.split(maxsplit=1)[0]
        if name1 and name2 and name1[-1] in '123' and name2[-1] in '123':
            return name1[:-1] == name2[:-1]
        return name1 == name2
    # Do not call .encode functions but use the unicode pointer inside the
    # python object directly, provided it is in 1-byte encoding, so we can
    # find the spaces and tabs easily.
    cdef char * header1_chars = <char *>PyUnicode_1BYTE_DATA(header1)
    cdef char * header2_chars = <char *>PyUnicode_1BYTE_DATA(header2)
    cdef size_t header1_length = <size_t>PyUnicode_GET_LENGTH(header1)
    return record_ids_match(header1_chars, header2_chars, header1_length)


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
