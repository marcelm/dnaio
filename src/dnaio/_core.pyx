# cython: language_level=3, emit_code_comments=False

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING
from libc.string cimport strncmp, memcmp, memcpy
cimport cython

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
                yield Sequence.__new__(Sequence, name, sequence, qualities)
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
    paired-end reads that have names ending in '/1' and '/2'. Also, the
    fastq-dump tool (used for converting SRA files to FASTQ) appends '.1', '.2'
    and sometimes '.3' to paired-end reads if option -I is used.
    """
    cdef:
        str name1, name2

    name1 = header1.split()[0]
    name2 = header2.split()[0]
    if name1 and name2 and name1[-1] in '123' and name2[-1] in '123':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2
