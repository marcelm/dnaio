# cython: language_level=3, emit_code_comments=False

from cpython.buffer cimport PyBUF_SIMPLE, PyObject_GetBuffer, PyBuffer_Release
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING, PyBytes_GET_SIZE, PyBytes_CheckExact
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
from cpython.unicode cimport PyUnicode_CheckExact, PyUnicode_GET_LENGTH, PyUnicode_DecodeASCII
from cpython.object cimport Py_TYPE, PyTypeObject
from cpython.ref cimport PyObject
from cpython.tuple cimport PyTuple_GET_ITEM
from libc.string cimport memcmp, memcpy, memchr, strcspn, strspn, memmove
cimport cython

cdef extern from "Python.h":
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)
    object PyUnicode_New(Py_ssize_t size, Py_UCS4 maxchar)

cdef extern from *:
    """
    #if defined(USE_SSE2)
      #include "ascii_check_sse2.h"
    #else
      #include "ascii_check.h"
    #endif
    """
    int string_is_ascii(char *string, size_t length)

cdef extern from "_conversions.h":
    const char NUCLEOTIDE_COMPLEMENTS[256]

from .exceptions import FastqFormatError
from ._util import shorten


def bytes_ascii_check(bytes string, Py_ssize_t length = -1):
    if length == -1:
        length = PyBytes_GET_SIZE(string)
    else:
        length = min(length, PyBytes_GET_SIZE(string))
    cdef bint ascii = string_is_ascii(PyBytes_AS_STRING(string), length)
    return ascii


def is_not_ascii_message(field, value):
    """
    Return an error message for a non-ASCII field encountered when initializing a SequenceRecord

    Arguments:
        field: Description of the field ("name", "sequence", "qualities" or similar)
            in which non-ASCII characters were found
        value: Unicode string that was intended to be assigned to the field
    """
    detail = ""
    try:
        value.encode("ascii")
    except UnicodeEncodeError as e:
        detail = (
            f", but found '{value[e.start:e.end]}' at index {e.start}"
        )
    return f"'{field}' in sequence file must be ASCII encoded{detail}"


cdef class SequenceRecord:
    """
    A named sequence with optional quality values.
    This typically represents a record from a FASTA or FASTQ file.
    The readers returned by `dnaio.open` yield objects of this type
    when mode is set to ``"r"``

    Attributes:
        name (str): The read header
        sequence (str): The nucleotide (or amino acid) sequence
        qualities (str): None if no quality values are available
            (such as when the record comes from a FASTA file).
            If quality values are available, this is a string
            that contains the Phred-scaled qualities encoded as
            ASCII(qual+33) (as in FASTQ).

    Raises:
        ValueError: One of the provide attributes is not ASCII or
            the lengths of sequence and qualities differ
    """
    cdef:
        object _name
        object _sequence
        object _qualities
        object _id
        object _comment

    def __init__(self, object name, object sequence, object qualities = None):
        if not PyUnicode_CheckExact(name):
            raise TypeError(f"name should be of type str, got {type(name)}")
        if not PyUnicode_IS_COMPACT_ASCII(name):
            raise ValueError(is_not_ascii_message("name", name))
        if not PyUnicode_CheckExact(sequence):
            raise TypeError(f"sequence should be of type str, got {type(sequence)}")
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError(is_not_ascii_message("sequence", sequence))
        if qualities is not None:
            if not PyUnicode_CheckExact(qualities):
                raise TypeError(f"qualities should be of type str, got {type(qualities)}")
            if not PyUnicode_IS_COMPACT_ASCII(qualities):
                raise ValueError(is_not_ascii_message("qualities", qualities))
            if len(qualities) != len(sequence):
                rname = shorten(name)
                raise ValueError("In read named {!r}: length of quality sequence "
                                 "({}) and length of read ({}) do not match".format(
                    rname, len(qualities), len(sequence)))
        self._name = name
        self._sequence = sequence
        self._qualities = qualities

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if not PyUnicode_CheckExact(name):
            raise TypeError(f"name must be of type str, got {type(name)}")
        if not PyUnicode_IS_COMPACT_ASCII(name):
            raise ValueError(is_not_ascii_message("name", name))
        self._name = name
        self._id = None
        self._comment = None

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        if not PyUnicode_CheckExact(sequence):
            raise TypeError(f"sequence must be of type str, got {type(sequence)}")
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError(is_not_ascii_message("sequence", sequence))
        self._sequence = sequence

    @property
    def qualities(self):
        return self._qualities

    @qualities.setter
    def qualities(self, qualities):
        if PyUnicode_CheckExact(qualities):
            if not PyUnicode_IS_COMPACT_ASCII(qualities):
                raise ValueError(is_not_ascii_message("qualities", qualities))
        elif qualities is None:
            pass
        else:
            raise TypeError(
                f"qualities must be of type str or None, "
                f"got {type(qualities)}."
            )
        self._qualities = qualities

    @property
    def id(self):
        cdef char *name
        cdef size_t name_length
        cdef size_t id_length
        # Not yet cached is None
        if self._id is None:
            name = <char *>PyUnicode_DATA(self._name)
            name_length = <size_t>PyUnicode_GET_LENGTH(self._name)
            id_length = strcspn(name, "\t ")
            if id_length == name_length:
                self._id = self._name
            else:
                self._id = PyUnicode_New(id_length, 127)
                memcpy(PyUnicode_DATA(self._id), name, id_length)
        return self._id

    @property
    def comment(self):
        cdef char *name
        cdef size_t name_length
        cdef size_t id_length
        cdef char *comment_start
        cdef size_t comment_length
        # Not yet cached is None
        if self._comment is None:
            name = <char *>PyUnicode_DATA(self._name)
            name_length = <size_t>PyUnicode_GET_LENGTH(self._name)
            id_length = strcspn(name, "\t ")
            if id_length == name_length:
                self._comment = ""
            else:
                comment_start = name + id_length + 1
                # Skip empty whitespace before comment
                comment_start = comment_start + strspn(comment_start, '\t ')
                comment_length = name_length - (comment_start - name)
                self._comment = PyUnicode_New(comment_length , 127)
                memcpy(PyUnicode_DATA(self._comment), comment_start, comment_length)
        # Cached but nothing is internally empty string, expose externally as None
        if PyUnicode_GET_LENGTH(self._comment) == 0:
            return None
        return self._comment

    def __getitem__(self, key):
        """
        Slice this SequenceRecord. If the qualities attribute is not None, it is
        sliced accordingly. The read name is copied unchanged.

        Returns:
            A new `SequenceRecord` object representing the sliced sequence.
        """
        return self.__class__(
            self._name,
            self._sequence[key],
            self._qualities[key] if self._qualities is not None else None,
        )

    def __repr__(self):
        qstr = ''
        if self._qualities is not None:
            qstr = ', qualities={!r}'.format(shorten(self._qualities))
        return '<SequenceRecord(name={!r}, sequence={!r}{})>'.format(
            shorten(self._name), shorten(self._sequence), qstr)

    def __len__(self):
        """
        Returns:
           The number of characters in the sequence
        """
        return len(self._sequence)

    def __richcmp__(self, SequenceRecord other, int op):
        if 2 <= op <= 3:
            eq = self._name == other._name and \
                self._sequence == other._sequence and \
                self._qualities == other._qualities
            if op == 2:
                return eq
            else:
                return not eq
        else:
            raise NotImplementedError()

    def __reduce__(self):
        return (SequenceRecord, (self._name, self._sequence, self._qualities))

    def qualities_as_bytes(self):
        """
        Return the qualities as a bytes object.

        This is a faster version of ``record.qualities.encode('ascii')``.
        """
        return self._qualities.encode('ascii')

    def fastq_bytes(self, bint two_headers=False):
        """
        Format this record in FASTQ format

        Arguments:
            two_headers (bool): If True, repeat the header (after the ``@``)
                on the third line (after the ``+``)

        Returns:
            A bytes object with the formatted record.
            This can be written directly to a file.
        """
        if self._qualities is None:
            raise ValueError("Cannot create a FASTQ record when qualities is not set.")

        cdef:
            char *name = <char *>PyUnicode_DATA(self._name)
            char *sequence = <char *>PyUnicode_DATA(self._sequence)
            char *qualities = <char *>PyUnicode_DATA(self._qualities)
            size_t name_length = <size_t>PyUnicode_GET_LENGTH(self._name)
            size_t sequence_length = <size_t>PyUnicode_GET_LENGTH(self._sequence)
            size_t qualities_length = <size_t>PyUnicode_GET_LENGTH(self._qualities)

        # Total size is name + sequence + qualities + 4 newlines + '+' and an
        # '@' to be put in front of the name.
        cdef Py_ssize_t total_size = name_length + sequence_length + qualities_length + 6

        if two_headers:
            # We need space for the name after the +.
            total_size += name_length

        # This is the canonical way to create an uninitialized bytestring of given size
        cdef bytes retval = PyBytes_FromStringAndSize(NULL, total_size)
        cdef char *retval_ptr = PyBytes_AS_STRING(retval)

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


    def fastq_bytes_two_headers(self):
        # Deprecated, use ``.fastq_bytes(two_headers=True)`` instead.
        return self.fastq_bytes(two_headers=True)

    def is_mate(self, SequenceRecord other):
        """
        Check whether this instance and another are part of the same read pair

        Checking is done by comparing IDs. The ID is the part of the name
        before the first whitespace. Any 1, 2 or 3 at the end of the IDs is
        excluded from the check as forward reads may have a 1 appended to their
        ID and reverse reads a 2 etc.

        Args:
            other (SequenceRecord): The object to compare to

        Returns:
            bool: Whether this and *other* are part of the same read pair.
        """
        cdef:
            char *header1_chars = <char *>PyUnicode_DATA(self._name)
            char *header2_chars = <char *>PyUnicode_DATA(other._name)
            size_t header2_length = <size_t>PyUnicode_GET_LENGTH(other._name)
            size_t id1_length = strcspn(header1_chars, ' \t')
            bint id1_ends_with_number = b'1' <= header1_chars[id1_length - 1] <= b'3'
        return record_ids_match(header1_chars, header2_chars, id1_length,
                                header2_length, id1_ends_with_number)

    def reverse_complement(self):
        """
        Return a reverse-complemented version of this record.

        - The name remains unchanged.
        - The sequence is reverse complemented.
        - If quality values exist, their order is reversed.
        """
        cdef:
            Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(self._sequence)
            object reversed_sequence_obj = PyUnicode_New(sequence_length, 127)
            object reversed_qualities_obj
            char *reversed_sequence = <char *>PyUnicode_DATA(reversed_sequence_obj)
            char *sequence = <char *>PyUnicode_DATA(self._sequence),
            char *reversed_qualities
            char *qualities
            Py_ssize_t cursor, reverse_cursor
            unsigned char nucleotide
            SequenceRecord seq_record
        reverse_cursor = sequence_length
        for cursor in range(sequence_length):
            reverse_cursor -= 1
            nucleotide = sequence[cursor]
            reversed_sequence[reverse_cursor] = NUCLEOTIDE_COMPLEMENTS[nucleotide]

        if self._qualities is not None:
            reverse_cursor = sequence_length
            reversed_qualities_obj = PyUnicode_New(sequence_length, 127)
            reversed_qualities = <char *>PyUnicode_DATA(reversed_qualities_obj)
            qualities = <char *>PyUnicode_DATA(self._qualities)
            for cursor in range(sequence_length):
                reverse_cursor -= 1
                reversed_qualities[reverse_cursor] = qualities[cursor]
        else:
            reversed_qualities_obj = None
        seq_record = SequenceRecord.__new__(SequenceRecord)
        seq_record._name = self._name
        seq_record._sequence = reversed_sequence_obj
        seq_record._qualities = reversed_qualities_obj
        return seq_record


def paired_fastq_heads(buf1, buf2, Py_ssize_t end1, Py_ssize_t end2):
    """
    Skip forward in the two buffers by multiples of four lines.

    Returns:
        A tuple (length1, length2) such that buf1[:length1] and
        buf2[:length2] contain the same number of lines (where the
        line number is divisible by four).
    """
    # Acquire buffers. Cython automatically checks for errors here.
    cdef Py_buffer data1_buffer
    cdef Py_buffer data2_buffer
    PyObject_GetBuffer(buf1, &data1_buffer, PyBUF_SIMPLE)
    PyObject_GetBuffer(buf2, &data2_buffer, PyBUF_SIMPLE)

    cdef:
        Py_ssize_t linebreaks = 0
        char *data1 = <char *>data1_buffer.buf
        char *data2 = <char *>data2_buffer.buf
        # The min() function ensures we do not read beyond the size of the buffer.
        char *data1_end = data1 + min(end1, data1_buffer.len)
        char *data2_end = data2 + min(end2, data2_buffer.len)
        char *pos1 = data1
        char *pos2 = data2
        char *record_start1 = data1
        char *record_start2 = data2

    while True:
        pos1 = <char *>memchr(pos1, b'\n', data1_end - pos1)
        if pos1 == NULL:
            break
        pos1 += 1
        pos2 = <char *>memchr(pos2, b'\n', data2_end - pos2)
        if pos2 == NULL:
            break
        pos2 += 1
        linebreaks += 1
        if linebreaks == 4:
            linebreaks = 0
            record_start1 = pos1
            record_start2 = pos2

    # Hit the end of the data block
    # This code will always be reached, so the buffers are always safely released.
    PyBuffer_Release(&data1_buffer)
    PyBuffer_Release(&data2_buffer)
    return record_start1 - data1, record_start2 - data2


cdef class FastqIter:
    """
    Parse a FASTQ file and yield SequenceRecord objects

    Arguments:
        file: a file-like object, opened in binary mode (it must have a readinto
            method)

        sequence_class: A custom class to use for the returned instances
            (instead of SequenceRecord)

        buffer_size: size of the initial buffer. This is automatically grown
            if a FASTQ record is encountered that does not fit.

    Yields:
        The *first value* that the generator yields is a boolean indicating whether
        the first record in the FASTQ has a repeated header (in the third row
        after the ``+``). Subsequent values are SequenceRecord objects (or whichever
        objects sequence_class returned if specified)
    """
    cdef:
        Py_ssize_t buffer_size
        char *buffer
        Py_ssize_t bytes_in_buffer
        type sequence_class
        bint use_custom_class
        bint extra_newline
        bint yielded_two_headers
        bint eof
        object file
        char *record_start
    cdef readonly Py_ssize_t number_of_records

    def __cinit__(self, file, sequence_class, Py_ssize_t buffer_size):
        self.buffer_size = buffer_size
        self.buffer = <char *>PyMem_Malloc(self.buffer_size)
        if self.buffer == NULL:
            raise MemoryError()
        self.bytes_in_buffer = 0
        self.sequence_class = sequence_class
        self.use_custom_class = sequence_class is not SequenceRecord
        self.number_of_records = 0
        self.extra_newline = False
        self.yielded_two_headers = False
        self.eof = False
        self.record_start = self.buffer
        self.file = file
        if buffer_size < 1:
            raise ValueError("Starting buffer size too small")

    def __dealloc__(self):
        PyMem_Free(self.buffer)

    cdef _read_into_buffer(self):
        # This function sets self.record_start at 0 and makes sure self.buffer
        # starts at the start of a FASTQ record. Any incomplete FASTQ remainder
        # of the already processed buffer is moved to the start of the buffer
        # and the rest of the buffer is filled up with bytes from the file.

        cdef char *tmp
        cdef Py_ssize_t remaining_bytes

        if self.record_start == self.buffer and self.bytes_in_buffer == self.buffer_size:
            # buffer too small, double it
            self.buffer_size *= 2
            tmp = <char *>PyMem_Realloc(self.buffer, self.buffer_size)
            if tmp == NULL:
                raise MemoryError()
            self.buffer = tmp
        else:
            # Move the incomplete record from the end of the buffer to the beginning.
            remaining_bytes = self.bytes_in_buffer - (self.record_start - self.buffer)
            # Memmove copies safely when dest and src overlap.
            memmove(self.buffer, self.record_start, remaining_bytes)
            self.bytes_in_buffer = remaining_bytes
        self.record_start = self.buffer

        cdef Py_ssize_t empty_bytes_in_buffer = self.buffer_size - self.bytes_in_buffer
        cdef object filechunk = self.file.read(empty_bytes_in_buffer)
        if not PyBytes_CheckExact(filechunk):
            raise TypeError("self.file is not a binary file reader.")
        cdef Py_ssize_t filechunk_size = PyBytes_GET_SIZE(filechunk)
        if filechunk_size > empty_bytes_in_buffer:
            raise ValueError(f"read() returned too much data: "
                             f"{empty_bytes_in_buffer} bytes requested, "
                             f"{filechunk_size} bytes returned.")
        memcpy(self.buffer + self.bytes_in_buffer, PyBytes_AS_STRING(filechunk), filechunk_size)
        # Strings are tested for ASCII as FASTQ should only contain ASCII characters.
        if not string_is_ascii(self.buffer + self.bytes_in_buffer,
                               filechunk_size):
            raise FastqFormatError(
                "Non-ASCII characters found in record.", None)
        self.bytes_in_buffer += filechunk_size

        if filechunk_size == 0:  # End of file
            if self.bytes_in_buffer == 0:  # EOF Reached. Stop iterating.
                self.eof = True
            elif not self.extra_newline and self.buffer[self.bytes_in_buffer - 1] != b'\n':
                # There is still data in the buffer and its last character is
                # not a newline: This is a file that is missing the final
                # newline. Append a newline and continue.
                self.buffer[self.bytes_in_buffer] = b'\n'
                self.bytes_in_buffer += 1
                self.extra_newline = True
            else:  # Incomplete FASTQ records are present.
                if self.extra_newline:
                    # Do not report the linefeed that was added by dnaio but
                    # was not present in the original input.
                    self.bytes_in_buffer -= 1
                record = PyUnicode_DecodeASCII(self.record_start, self.bytes_in_buffer, NULL)
                lines = record.count('\n')
                raise FastqFormatError(
                    'Premature end of file encountered. The incomplete final record was: '
                    '{!r}'.format(shorten(record, 500)),
                    line=self.number_of_records * 4 + lines)

    def __iter__(self):
        return self

    def __next__(self):
        cdef:
            object ret_val
            SequenceRecord seq_record
            char *name_start
            char *name_end
            char *sequence_start
            char *sequence_end
            char *second_header_start
            char *second_header_end
            char *qualities_start
            char *qualities_end
            char *buffer_end
            size_t remaining_bytes
            Py_ssize_t name_length, sequence_length, second_header_length, qualities_length
        # Repeatedly attempt to parse the buffer until we have found a full record.
        # If an attempt fails, we read more data before retrying.
        while True:
            buffer_end = self.buffer + self.bytes_in_buffer
            if self.eof:
                raise StopIteration()
            ### Check for a complete record (i.e 4 newlines are present)
            # Use libc memchr, this optimizes looking for characters by
            # using 64-bit integers. See:
            # https://sourceware.org/git/?p=glibc.git;a=blob_plain;f=string/memchr.c;hb=HEAD
            # void *memchr(const void *str, int c, size_t n)
            name_end = <char *>memchr(self.record_start, b'\n', <size_t>(buffer_end - self.record_start))
            if name_end == NULL:
                self._read_into_buffer()
                continue
            # self.bytes_in_buffer - sequence_start is always nonnegative:
            # - name_end is at most self.bytes_in_buffer - 1
            # - thus sequence_start is at most self.bytes_in_buffer
            sequence_start = name_end + 1
            sequence_end = <char *>memchr(sequence_start, b'\n', <size_t>(buffer_end - sequence_start))
            if sequence_end == NULL:
                self._read_into_buffer()
                continue
            second_header_start = sequence_end + 1
            remaining_bytes = (buffer_end - second_header_start)
            # Usually there is no second header, so we skip the memchr call.
            if remaining_bytes > 2 and second_header_start[0] == b'+' and second_header_start[1] == b'\n':
                second_header_end = second_header_start + 1
            else:
                second_header_end = <char *>memchr(second_header_start, b'\n', <size_t>(remaining_bytes))
                if second_header_end == NULL:
                    self._read_into_buffer()
                    continue
            qualities_start = second_header_end + 1
            qualities_end = <char *>memchr(qualities_start, b'\n', <size_t>(buffer_end - qualities_start))
            if qualities_end == NULL:
                self._read_into_buffer()
                continue

            if self.record_start[0] != b'@':
                raise FastqFormatError("Line expected to "
                    "start with '@', but found {!r}".format(chr(self.record_start[0])),
                    line=self.number_of_records * 4)
            if second_header_start[0] != b'+':
                raise FastqFormatError("Line expected to "
                    "start with '+', but found {!r}".format(chr(second_header_start[0])),
                    line=self.number_of_records * 4 + 2)

            name_start = self.record_start + 1  # Skip @
            second_header_start += 1  # Skip +
            name_length = name_end - name_start
            sequence_length = sequence_end - sequence_start
            second_header_length = second_header_end - second_header_start
            qualities_length = qualities_end - qualities_start

            # Check for \r\n line-endings and compensate
            if (name_end - 1)[0] == b'\r':
                name_length -= 1
            if (sequence_end - 1)[0] == b'\r':
                sequence_length -= 1
            if (second_header_end - 1)[0] == b'\r':
                second_header_length -= 1
            if (qualities_end - 1)[0] == b'\r':
                qualities_length -= 1

            if second_header_length:  # should be 0 when only + is present
                if (name_length != second_header_length or
                        memcmp(second_header_start, name_start, second_header_length) != 0):
                    raise FastqFormatError(
                        "Sequence descriptions don't match ('{}' != '{}').\n"
                        "The second sequence description must be either "
                        "empty or equal to the first description.".format(
                            PyUnicode_DecodeASCII(name_start, name_length, NULL),
                            PyUnicode_DecodeASCII(second_header_start, second_header_length, NULL)),
                        line=self.number_of_records * 4 + 2)

            if qualities_length != sequence_length:
                raise FastqFormatError(
                    "Length of sequence and qualities differ", line=self.number_of_records * 4 + 3)

            if self.number_of_records == 0 and not self.yielded_two_headers:
                self.yielded_two_headers = True
                return bool(second_header_length)  # first yielded value is special

            # Constructing objects with PyUnicode_New and memcpy bypasses some of
            # the checks otherwise done when using PyUnicode_DecodeLatin1 or similar
            name = PyUnicode_New(name_length, 127)
            sequence = PyUnicode_New(sequence_length, 127)
            qualities = PyUnicode_New(qualities_length, 127)
            if <PyObject*>name == NULL or <PyObject*>sequence == NULL or <PyObject*>qualities == NULL:
                raise MemoryError()
            memcpy(PyUnicode_DATA(name), name_start, name_length)
            memcpy(PyUnicode_DATA(sequence), sequence_start, sequence_length)
            memcpy(PyUnicode_DATA(qualities), qualities_start, qualities_length)

            if self.use_custom_class:
                ret_val = self.sequence_class(name, sequence, qualities)
            else:
                seq_record = SequenceRecord.__new__(SequenceRecord)
                seq_record._name = name
                seq_record._sequence = sequence
                seq_record._qualities = qualities
                ret_val = seq_record
            # Advance record to next position
            self.number_of_records += 1
            self.record_start = qualities_end + 1
            return ret_val


def record_names_match(header1: str, header2: str):
    """
    Check whether the sequence record ids id1 and id2 are compatible, ignoring a
    suffix of '1', '2' or '3'. This exception allows to check some old
    paired-end reads that have IDs ending in '/1' and '/2'. Also, the
    fastq-dump tool (used for converting SRA files to FASTQ) appends '.1', '.2'
    and sometimes '.3' to paired-end reads if option -I is used.

    Deprecated, use `SequenceRecord.is_mate` instead
    """
    cdef:
        char *header1_chars = NULL
        char *header2_chars = NULL
        size_t header1_length
    if PyUnicode_CheckExact(header1):
        if PyUnicode_IS_COMPACT_ASCII(header1):
            header1_chars = <char *>PyUnicode_DATA(header1)
        else:
            raise ValueError("header1 must be a valid ASCII-string.")
    else:
        raise TypeError(f"Header 1 is the wrong type. Expected bytes or string, "
                        f"got: {type(header1)}")

    if PyUnicode_CheckExact(header2):
        if PyUnicode_IS_COMPACT_ASCII(header2):
            header2_chars = <char *>PyUnicode_DATA(header2)
            header2_length = <size_t> PyUnicode_GET_LENGTH(header2)
        else:
            raise ValueError("header2 must be a valid ASCII-string.")
    else:
        raise TypeError(f"Header 2 is the wrong type. Expected bytes or string, "
                        f"got: {type(header2)}")
    cdef size_t id1_length = strcspn(header1_chars, ' \t')
    cdef bint id1_ends_with_number =  b'1' <= header1_chars[id1_length - 1] <= b'3'
    return record_ids_match(header1_chars, header2_chars, id1_length,
                            header2_length, id1_ends_with_number)


cdef inline bint record_ids_match(char *header1,
                                  char *header2,
                                  size_t id1_length,
                                  size_t header2_length,
                                  bint id1_ends_with_number):
    """
    Check whether the ASCII-encoded IDs match.

    header1, header2 pointers to the ASCII-encoded headers
    id1_length, the length of header1 before the first whitespace
    header2_length, the full length of header2.
    id1_ends_with_number, whether id1 ends with a 1,2 or 3.
    """

    if header2_length < id1_length:
        return False

    cdef char end = header2[id1_length]
    if end != b'\000' and end != b' ' and end != b'\t':
        return False

    # Check if the IDs end with 1, 2 or 3. This is the read pair number
    # which should not be included in the comparison.
    cdef bint id2_ends_with_number = b'1' <= header2[id1_length - 1] <= b'3'
    if id1_ends_with_number and id2_ends_with_number:
        id1_length -= 1

    # Compare the strings up to the ID end position.
    return memcmp(<void *>header1, <void *>header2, id1_length) == 0


def records_are_mates(*args) -> bool:
    """
    Check if the provided `SequenceRecord` objects are all mates of each other by
    comparing their record IDs.
    Accepts two or more `SequenceRecord` objects.

    This is the same as `SequenceRecord.is_mate` in the case of only two records,
    but allows for for cases where information is split into three records or more
    (such as UMI, R1, R2 or index, R1, R2).

    If there are only two records to check, prefer `SequenceRecord.is_mate`.

    Example usage::

        for records in zip(*all_my_fastq_readers):
            if not records_are_mates(*records):
                raise MateError(f"IDs do not match for {records}")

    Args:
        *args: two or more `~dnaio.SequenceRecord` objects

    Returns: True or False
    """
    cdef Py_ssize_t args_length = len(args)
    if args_length < 2:
        raise TypeError("records_are_mates requires at least two arguments")

    cdef SequenceRecord first = <SequenceRecord>PyTuple_GET_ITEM(args, 0)
    if Py_TYPE(first) != <PyTypeObject *> SequenceRecord:
        raise TypeError(f"{first:r} is not a SequenceRecord object")

    cdef:
        object first_name_obj = first._name
        char *first_name = <char *>PyUnicode_DATA(first_name_obj)
        Py_ssize_t first_name_length = PyUnicode_GET_LENGTH(first_name_obj)
        Py_ssize_t id_length = strcspn(first_name, b' \t')
        bint id_ends_with_number = b'1' <= first_name[id_length - 1] <= b'3'
        SequenceRecord other
        object other_name_obj
        char *other_name
        Py_ssize_t other_name_length
        bint other_id_ends_with_number
        char end_char
        bint are_mates = True
        Py_ssize_t i

    for i in range(1, args_length):
        other = <SequenceRecord>PyTuple_GET_ITEM(args, i)
        if Py_TYPE(other) != <PyTypeObject *>SequenceRecord:
            raise TypeError(f"{other:r} is not a SequenceRecord object")
        other_name_obj = other._name
        other_name = <char *>PyUnicode_DATA(other_name_obj)
        other_name_length = PyUnicode_GET_LENGTH(other_name_obj)
        # If a match is false, are_mates will stay false regardless of any true checks afterward.
        are_mates &= record_ids_match(first_name, other_name, id_length,
                                      other_name_length, id_ends_with_number)
    return are_mates
