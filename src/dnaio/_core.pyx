# cython: language_level=3, emit_code_comments=False
from cpython.buffer cimport PyBUF_SIMPLE, PyObject_GetBuffer, PyBuffer_Release
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING, PyBytes_GET_SIZE, PyBytes_CheckExact, _PyBytes_Resize
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
from cpython.number cimport PyNumber_AsSsize_t
from cpython.unicode cimport PyUnicode_CheckExact, PyUnicode_GET_LENGTH, PyUnicode_DecodeASCII
from cpython.object cimport Py_TYPE, PyTypeObject
from cpython.pyport cimport PY_SSIZE_T_MAX
from cpython.ref cimport PyObject
from cpython.slice cimport PySlice_Check, PySlice_GetIndicesEx
from cpython.tuple cimport PyTuple_GET_ITEM
from libc.string cimport memcmp, memcpy, memchr, strcspn, strspn, memmove, strlen
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int8_t, int16_t, int32_t, \
    UINT8_MAX, UINT16_MAX, UINT32_MAX, INT8_MIN, INT8_MAX, INT16_MIN, INT16_MAX, \
    INT32_MIN, INT32_MAX

cimport cython

cdef extern from "Python.h":
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)
    object PyUnicode_New(Py_ssize_t size, Py_UCS4 maxchar)

cdef extern from "ascii_check.h":
    int string_is_ascii(char *string, size_t length)

cdef extern from "_conversions.h":
    const char NUCLEOTIDE_COMPLEMENTS[256]

cdef extern from "bam.h":
    void decode_bam_sequence(void *dest, void *encoded_sequence, size_t length)
    void decode_bam_qualities(uint8_t *dest, uint8_t *encoded_qualities, size_t length)

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

cdef inline Py_ssize_t get_tag_int_value(uint8_t *tag_ptr):
    """
    Retrieve a integer value from a known integer tag. Return PY_SSIZE_T_MAX
    when an error occurs.

    There is a tendency for BAM programs to choose the smallest possible
    representation for an integer. So even for known integer tags it is unsure
    what the type is going to be.
    """
    cdef uint8_t tag_type = tag_ptr[2]
    cdef uint8_t *value_ptr = tag_ptr + 3
    if tag_type == b"c":
        return (<int8_t *> value_ptr)[0]
    elif tag_type == b"C":
        return (<uint8_t *>value_ptr)[0]
    elif tag_type == b"s":
        return (<int16_t *>value_ptr)[0]
    elif tag_type == b"S":
        return (<uint16_t *>value_ptr)[0]
    elif tag_type == b"i":
        return (<int32_t *>value_ptr)[0]
    elif tag_type == b"I":
        return (<uint32_t *>value_ptr)[0]
    else:
        return PY_SSIZE_T_MAX

cdef inline Py_ssize_t store_tag_int_value(char *tag, Py_ssize_t value, uint8_t *dest):
    """
    Writes both the tag type and the value at dest.
    Returns the amount of written bytes or -1 for error.
    """
    cdef void *value_ptr = <void *>(dest + 3)
    if value < INT32_MIN or value > UINT32_MAX:
        return -1
    memcpy(dest, tag, 2)
    if INT8_MIN <= value <= INT8_MAX:
       dest[2] = b'c'
       (<int8_t *>value_ptr)[0] = <int8_t>value
       return 4
    elif INT8_MAX < value <= UINT8_MAX:
       dest[2] = b'C'
       (<uint8_t *>value_ptr)[0] = <uint8_t>value
       return 4
    if INT16_MIN <= value <= INT16_MAX:
       dest[2] = b's'
       (<int16_t *>value_ptr)[0] = <int16_t>value
       return 5
    elif INT16_MAX < value <= UINT16_MAX:
       dest[2] = b'S'
       (<uint16_t *>value_ptr)[0] = <uint16_t>value
       return 5
    if INT32_MIN <= value <= INT32_MAX:
       dest[2] = b'i'
       (<int32_t *>value_ptr)[0] = <int32_t>value
       return 7
    elif INT32_MAX < value <= UINT32_MAX:
       dest[2] = b'I'
       (<uint32_t *>value_ptr)[0] = <uint32_t>value
       return 7
    return -1


cdef inline Py_ssize_t get_array_item_size(uint8_t array_type):
    if array_type == b'c' or array_type == b'C':
        return 1
    elif array_type == b's' or array_type == b'S':
        return 2
    elif array_type == b'i' or array_type == b'I' or array_type == b'f':
        return 4
    else:
        return -1


cdef inline Py_ssize_t get_tag_length(uint8_t *tag, Py_ssize_t max_length):
    cdef:
        uint8_t tag_type = tag[2]
        uint8_t *tag_end = NULL
        Py_ssize_t array_length
        Py_ssize_t item_size
        Py_ssize_t tag_length = -1
    if tag_type == b"c" or tag_type == b"C" or tag_type == b"A":
        tag_length = 4
    elif tag_type == b"s" or tag_type == b"S":
        tag_length = 5
    elif tag_type == b"i" or tag_type == b"I" or tag_type == b"F":
        tag_length = 7
    elif tag_type == b"H" or b"Z":
        tag_end = <uint8_t *>memchr(<void *> tag + 3, 0, max_length)
        if tag_end == NULL:
            return -1
    elif tag_type == b"B":
        item_size = get_array_item_size(tag[3])
        if item_size == -1:
            return -1
        array_length = (<uint32_t *>(tag + 4))[0]
        tag_length = 8 + (item_size * array_length)
    else:
        tag_length = -1
    if tag_length > max_length:
        return -1
    return tag_length


cdef unpack_index_key(
    object key, Py_ssize_t original_size, Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step, Py_ssize_t *slice_length):
    """An index can be a slice or a number. Unpack both."""
    if PySlice_Check(key):
        PySlice_GetIndicesEx(
            key, original_size, start, stop, step, slice_length)
        return
    cdef:
        Py_ssize_t slice_step = 1
        Py_ssize_t slice_start = PyNumber_AsSsize_t(key, OverflowError)
        Py_ssize_t slice_stop
        Py_ssize_t length = 1
    if slice_start < 0:
        slice_start = original_size - slice_start
        if slice_start < 0:
            slice_start = 0
    if slice_start >= original_size:
        slice_start = original_size
        slice_stop = original_size
        length = 0
    else:
        slice_stop = slice_start + 1
        length = 1
    start[0] = slice_start
    stop[0] = slice_stop
    step[0] = slice_step
    slice_length[0] = length


cdef Py_ssize_t trim_move_table_tag(
    uint8_t *mv_tag,
    uint8_t *dest,
    Py_ssize_t ts,
    Py_ssize_t ns,
    Py_ssize_t start,
    Py_ssize_t stop):
    """
    Trims the move table in the mv_tag and recalculates the ts and ns values.
    the new mv, ts and ns tags are stored in dest.
    return the amount of bytes written to dest
    """
    if memcmp(mv_tag, b"mvBc", 4) != 0:
        return -1
    cdef:
        size_t raw_length = (<uint32_t *>(mv_tag + 4))[0]
        size_t mv_length = raw_length - 1
        # First entry in the mv_table is the stride.
        size_t stride = mv_tag[8]
        # Actual move table starts after that at a +1 offset.
        int8_t *mv_table = <int8_t *>(mv_tag + 9)
        int8_t *cursor = mv_table
        int8_t *mv_table_end = mv_table + mv_length
        Py_ssize_t sequence_pos = 0
        Py_ssize_t trimmed_positions = 0
        int8_t mv
    # Find the beginning
    while cursor < mv_table_end:
        mv = cursor[0]
        if mv == 1:
            sequence_pos += 1
        # 1,0,0,0,1,0,0,0,0,1,0,0,1,0
        # 0 starts here
        #         1 starts here
        #                   2 starts here
        if sequence_pos > start:
            sequence_pos -= 1
            break
        trimmed_positions += 1
        cursor += 1
    cdef int8_t *new_mv_table_start = cursor
    while cursor < mv_table_end and sequence_pos < stop:
        mv = cursor[0]
        if mv == 1:
            sequence_pos += 1
        cursor += 1
    cdef int8_t *new_mv_table_end = cursor
    cdef Py_ssize_t new_mv_table_size = new_mv_table_end - new_mv_table_start
    if new_mv_table_size == 0:
        ts = -1
        ns = -1
    else:
        if ts >= 0:
            ts += (trimmed_positions * stride)
        if ns >= 0:
            # Treat ts as 0 if missing (-1)
            ns = new_mv_table_size * stride + max(ts, 0)
    cdef uint8_t *dest_ptr = dest
    memcpy(dest, b"mvBc", 4)
    dest_ptr += 4
    (<uint32_t *>dest_ptr)[0] = new_mv_table_size + 1
    dest_ptr += 4
    dest_ptr[0] = stride
    dest_ptr += 1
    memcpy(dest_ptr, new_mv_table_start, new_mv_table_size)
    dest_ptr += new_mv_table_size
    cdef Py_ssize_t ns_tag_size = store_tag_int_value("ns", ns, dest_ptr)
    if ns_tag_size >= 0:
        dest_ptr += ns_tag_size
    cdef Py_ssize_t ts_tag_size = store_tag_int_value("ts", ts, dest_ptr)
    if ts_tag_size >= 0:
        dest_ptr += ts_tag_size
    return dest_ptr - dest

cdef Py_ssize_t trim_mm_ml_tag(
    uint8_t *MM,
    uint8_t *ML,
    uint8_t *dest,
    Py_ssize_t start,
    Py_ssize_t stop
):
    return 0


cdef inline object slice_tags(bytes tags, Py_ssize_t original_size, slice_obj):
    cdef:
        uint8_t *tag_start = <uint8_t *>PyBytes_AS_STRING(tags)
        uint8_t *tag = tag_start
        Py_ssize_t tags_length = PyBytes_GET_SIZE(tags)
        uint8_t *tags_end = tag_start + tags_length
        Py_ssize_t start = 0
        Py_ssize_t stop = original_size
        Py_ssize_t step = 1
        Py_ssize_t slice_length = original_size

    unpack_index_key(slice_obj, original_size, &start, &stop, &step, &slice_length)
    if start == 0 and stop == original_size:
        # No need to do work
        return tags
    cdef:
        object sliced_tags = PyBytes_FromStringAndSize(NULL, tags_length)
        PyObject *sliced_tags_ptr = <PyObject *>sliced_tags
        uint8_t *destination = <uint8_t *>PyBytes_AS_STRING(sliced_tags)
        uint8_t *dest_ptr = destination
        Py_ssize_t ns = -1
        Py_ssize_t ts = -1
        uint8_t *MM = NULL
        uint8_t *ML = NULL
        uint8_t *mv = NULL

    while tag < tags_end:
        # Do not directly write ns, ts and MN as these are invalidated by cutting.
        tag_length  = get_tag_length(tag, tags_end - tag)
        if tag_length == -1:
            raise ValueError(f"Invalid tag starting from {tags[tag - tag_start:]}")
        if memcmp(tag, b"ns", 2) == 0:
            ns = get_tag_int_value(tag)
            if ns == PY_SSIZE_T_MAX:
                ns = -1  # Non-dorado ns tag. Just ignore.
            tag += tag_length
            continue
        if memcmp(tag, b"ts", 2) == 0:
            ts = get_tag_int_value(tag)
            if ts == PY_SSIZE_T_MAX:
                ts = -1  # Non-dorado ts tag. Just ignore.
            tag += tag_length
            continue
        if memcmp(tag, b"mv", 2) == 0:  # Signal to move table
            mv = tag
            tag += tag_length
            continue
        if memcmp(tag, b"ML", 2) == 0:
            ML = tag
            tag += tag_length
            continue
        if memcmp(tag, b"MM", 2) == 0:
            MM = tag
            tag += tag_length
            continue
        if memcmp(tag, b"MN", 2) == 0:
            # Value is invalidated by slicing.
            tag += tag_length
            continue
        if (memcmp(tag, b"du", 2) == 0):
            # Skip duration tag, it only makes sense for original length.
            tag += tag_length
            continue
        # In other cases just copy the tag
        memcpy(dest_ptr, tag, tag_length)
        tag += tag_length
        dest_ptr += tag_length

    cdef Py_ssize_t trim_move_table_size = 0
    cdef Py_ssize_t trim_mm_ml_size = 0
    if step == 1:
        # If step != 1 the appropriate cutting is not implemented. So the tags
        # are simply removed instead of adapted.
        trim_move_table_size = trim_move_table_tag(mv, dest_ptr, ts, ns, start, stop)
        if trim_move_table_size >= 0:
            dest_ptr += trim_move_table_size
        trim_mm_ml_size = trim_mm_ml_tag(MM, ML, dest_ptr, start, stop)
        if trim_mm_ml_size >= 0:
            dest_ptr += trim_mm_ml_size
    cdef Py_ssize_t destination_size = dest_ptr - destination
    if destination_size != tags_length:
        _PyBytes_Resize(&sliced_tags_ptr, destination_size)
    return <object>sliced_tags_ptr


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
        readonly object _tags

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
        """
        The header part before any whitespace. This is the unique identifier
        for the sequence.
        """
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
        """
        The header part after the first whitespace. This is usually used
        to store metadata. It may be empty in which case the attribute is None.
        """
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
        # Empty comment is returned as None. This is not stored internally as
        # None, otherwise the above code would run every time the attribute
        # was accessed.
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
            slice_tags(self._tags, PyUnicode_GET_LENGTH(self._sequence), key) if self._tags is not None else None,
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
            if remaining_bytes > 2 and memcmp(second_header_start, b"+\n", 2) == 0:
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

cdef struct BamRecordHeader:
    uint32_t block_size
    int32_t reference_id
    int32_t pos
    uint8_t l_read_name
    uint8_t mapq
    uint16_t bin
    uint16_t n_cigar_op
    uint16_t flag
    uint32_t l_seq
    int32_t next_ref_id
    int32_t next_pos
    int32_t tlen

cdef class BamIter:
    cdef:
        uint8_t *record_start
        uint8_t *buffer_end
        size_t read_in_size
        uint8_t *read_in_buffer
        size_t read_in_buffer_size
        object file
        readonly object header
        readonly Py_ssize_t number_of_records

    def __dealloc__(self):
        PyMem_Free(self.read_in_buffer)

    def __cinit__(self, fileobj, read_in_size = 48 * 1024):
        if read_in_size < 4:
            raise ValueError(f"read_in_size must be at least 4 got "
                              f"{read_in_size}")

        # Skip ahead and save the BAM header for later inspection
        magic_and_header_size = fileobj.read(8)
        if not isinstance(magic_and_header_size, bytes):
            raise TypeError(f"fileobj {fileobj} is not a binary IO type, "
                            f"got {type(fileobj)}")
        if len(magic_and_header_size) < 8:
            raise EOFError("Truncated BAM file")
        if magic_and_header_size[:4] != b"BAM\1":
            raise ValueError(
                f"fileobj: {fileobj}, is not a BAM file. No BAM magic, instead "
                f"found {magic_and_header_size[:4]}")
        l_text = int.from_bytes(magic_and_header_size[4:], "little", signed=False)
        header =  fileobj.read(l_text)
        if len(header) < l_text:
            raise EOFError("Truncated BAM file")
        n_ref_obj = fileobj.read(4)
        if len(n_ref_obj) < 4:
            raise EOFError("Truncated BAM file")
        n_ref = int.from_bytes(n_ref_obj, "little", signed=False)
        for i in range(n_ref):
            l_name_obj = fileobj.read(4)
            if len(l_name_obj) < 4:
                raise EOFError("Truncated BAM file")
            l_name = int.from_bytes(l_name_obj, "little", signed=False)
            reference_chunk_size = l_name + 4  # Include name and uint32_t of size
            reference_chunk = fileobj.read(reference_chunk_size)
            if len(reference_chunk) < reference_chunk_size:
                raise EOFError("Truncated BAM file")
        # Fileobj is now skipped ahead and at the start of the BAM records

        self.header = header
        self.read_in_size = read_in_size
        self.file = fileobj
        self.read_in_buffer = NULL
        self.read_in_buffer_size = 0
        self.record_start = self.read_in_buffer
        self.buffer_end = self.record_start

    def __iter__(self):
        return self

    cdef _read_into_buffer(self):
        cdef size_t read_in_size
        cdef size_t leftover_size = self.buffer_end - self.record_start
        cdef uint32_t block_size
        memmove(self.read_in_buffer, self.record_start, leftover_size)
        self.record_start = self.read_in_buffer
        self.buffer_end = self.record_start + leftover_size
        if leftover_size >= 4:
            # Immediately check how much data is at least required
            block_size = (<uint32_t *>self.record_start)[0]
            read_in_size = max(block_size, self.read_in_size)
        else:
            read_in_size = self.read_in_size - leftover_size
        new_bytes = self.file.read(read_in_size)
        cdef size_t new_bytes_size = PyBytes_GET_SIZE(new_bytes)
        cdef uint8_t *new_bytes_buf = <uint8_t *>PyBytes_AS_STRING(new_bytes)
        cdef size_t new_buffer_size = leftover_size + new_bytes_size
        if new_buffer_size == 0:
            # File completely read
            raise StopIteration()
        elif new_bytes_size == 0:
            raise EOFError("Incomplete record at the end of file")
        cdef uint8_t *tmp
        if new_buffer_size > self.read_in_buffer_size:
            tmp = <uint8_t *>PyMem_Realloc(self.read_in_buffer, new_buffer_size)
            if tmp == NULL:
                raise MemoryError()
            self.read_in_buffer = tmp
            self.read_in_buffer_size = new_buffer_size
        memcpy(self.read_in_buffer + leftover_size, new_bytes_buf, new_bytes_size)
        self.record_start = self.read_in_buffer
        self.buffer_end = self.record_start + new_buffer_size

    def __next__(self):
        cdef:
            SequenceRecord seq_record
            uint8_t *record_start
            uint8_t *buffer_end
            uint32_t record_size
            uint8_t *record_end
            cdef BamRecordHeader header
            cdef uint8_t *bam_name_start
            uint32_t name_length
            uint8_t *bam_seq_start
            uint32_t seq_length
            uint8_t *bam_qual_start
            uint8_t *bam_tags_start
            size_t bam_tags_length
            uint32_t encoded_seq_length

        while True:
            record_start = self.record_start
            buffer_end = self.buffer_end
            if record_start + 4 >= buffer_end:
                self._read_into_buffer()
                continue
            record_size = (<uint32_t *>record_start)[0]
            record_end = record_start + 4 + record_size
            if record_end > buffer_end:
                self._read_into_buffer()
                continue
            header = (<BamRecordHeader *>record_start)[0]
            if header.flag != 4:
                raise NotImplementedError(
                    "The BAM parser has been implemented with unmapped single "
                    "reads in mind to support ONT sequencing input. Mapped "
                    "BAM files or files with multiple reads are not supported. "
                    "Please use samtools fastq to make the appropriate "
                    "conversion to FASTQ format."
                )
            bam_name_start = record_start + sizeof(BamRecordHeader)
            name_length = header.l_read_name
            bam_seq_start = bam_name_start + name_length + header.n_cigar_op * sizeof(uint32_t)
            name_length -= 1  # Do not include the null byte
            seq_length = header.l_seq
            encoded_seq_length = (seq_length + 1) // 2
            bam_qual_start = bam_seq_start + encoded_seq_length
            bam_tags_start = bam_qual_start + seq_length
            bam_tags_length = record_end - bam_tags_start
            name = PyUnicode_New(name_length, 127)
            sequence = PyUnicode_New(seq_length, 127)
            qualities = PyUnicode_New(seq_length, 127)
            tags = PyBytes_FromStringAndSize(<char *>bam_tags_start, bam_tags_length)
            memcpy(<uint8_t *>PyUnicode_DATA(name), bam_name_start, name_length)
            decode_bam_sequence(<uint8_t *>PyUnicode_DATA(sequence), bam_seq_start, seq_length)
            if seq_length and bam_qual_start[0] == 0xff:
                # 0xff means qualities are missing.
                qualities = None
            else:
                qualities = PyUnicode_New(seq_length, 127)
                decode_bam_qualities(<uint8_t *>PyUnicode_DATA(qualities), bam_qual_start, seq_length)
            seq_record = SequenceRecord.__new__(SequenceRecord)
            seq_record._name = name
            seq_record._sequence = sequence
            seq_record._qualities = qualities
            seq_record._tags = tags
            self.number_of_records += 1
            self.record_start = record_end
            return seq_record

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
