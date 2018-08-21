"""
Classes for reading FASTA and FASTQ files
"""
import io
from xopen import xopen
from ._core import fastq_iter as _fastq_iter, Sequence, paired_fastq_heads as _paired_fastq_heads
from ._util import shorten as _shorten
from .exceptions import FileFormatError, UnknownFileFormat


class BinaryFileReader:
    """
    A mixin for readers that ensures that a file or a path can be passed in to the constructor.
    """
    _close_on_exit = False
    paired = False
    mode = 'rb'

    def __init__(self, file, _close_file=None):
        """
        The file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).
        """
        if isinstance(file, str):
            file = xopen(file, self.mode)
            self._close_on_exit = True
        elif _close_file:
            self._close_on_exit = True
        self._file = file

    def close(self):
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        if self._file is None:
            raise ValueError("I/O operation on closed BinaryFileReader")
        return self

    def __exit__(self, *args):
        self.close()


class FastaReader(BinaryFileReader):
    """
    Reader for FASTA files.
    """

    def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence, _close_file=None):
        """
        file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).

        keep_linebreaks -- whether to keep newline characters in the sequence
        """
        super().__init__(file, _close_file=_close_file)
        self.sequence_class = sequence_class
        self.delivers_qualities = False
        self._delimiter = '\n' if keep_linebreaks else ''

    def __iter__(self):
        """
        Read next entry from the file (single entry at a time).
        """
        name = None
        seq = []
        f = io.TextIOWrapper(self._file)
        for i, line in enumerate(f):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == '>':
                if name is not None:
                    yield self.sequence_class(name, self._delimiter.join(seq), None)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FileFormatError("At line {0}: Expected '>' at beginning of "
                    "FASTA record, but got {1!r}.".format(i + 1, _shorten(line)))

        if name is not None:
            yield self.sequence_class(name, self._delimiter.join(seq), None)
        # Prevent TextIOWrapper from closing the underlying file
        f.detach()


class FastqReader(BinaryFileReader):
    """
    Reader for FASTQ files. Does not support multi-line FASTQ files.
    """

    def __init__(self, file, sequence_class=Sequence, buffer_size=1048576, _close_file=None):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        """
        super().__init__(file, _close_file=_close_file)
        self.sequence_class = sequence_class
        self.delivers_qualities = True
        self.buffer_size = buffer_size

    def __iter__(self):
        return _fastq_iter(self._file, self.sequence_class, self.buffer_size)


def _fasta_head(buf, end):
    """
    Search for the end of the last complete FASTA record within buf[:end]

    Return an integer length such that buf[:length] contains the highest
    possible number of complete FASTA records.
    """
    pos = buf.rfind(b'\n>', 0, end)
    if pos != -1:
        return pos + 1
    if buf[0:1] == b'>':
        return 0
    raise FileFormatError('FASTA does not start with ">"')


def _fastq_head(buf, end=None):
    """
    Search for the end of the last complete *two* FASTQ records in buf[:end].

    Two FASTQ records are required to ensure that read pairs in interleaved
    paired-end data are not split.
    """
    linebreaks = buf.count(b'\n', 0, end)
    right = end
    for _ in range(linebreaks % 8 + 1):
        right = buf.rfind(b'\n', 0, right)
    # Note that this works even if linebreaks == 0:
    # rfind() returns -1 and adding 1 gives index 0,
    # which is correct.
    return right + 1


def read_chunks_from_file(f, buffer_size=4*1024**2):
    """
    Read a chunk of complete FASTA or FASTQ records from a file.
    The size of a chunk is at most buffer_size.
    f needs to be a file opened in binary mode.

    The yielded memoryview objects become invalid on the next iteration.
    """
    # This buffer is re-used in each iteration.
    buf = bytearray(buffer_size)

    # Read one byte to determine file format.
    # If there is a comment char, we assume FASTA!
    start = f.readinto(memoryview(buf)[0:1])
    if start == 1 and buf[0:1] == b'@':
        head = _fastq_head
    elif start == 1 and (buf[0:1] == b'#' or buf[0:1] == b'>'):
        head = _fasta_head
    else:
        raise UnknownFileFormat('Input file format unknown')

    # Layout of buf
    #
    # |-- complete records --|
    # +---+------------------+---------+-------+
    # |   |                  |         |       |
    # +---+------------------+---------+-------+
    # ^   ^                   ^         ^       ^
    # 0   start               end       bufend  len(buf)
    #
    # buf[0:start] is the 'leftover' data that could not be processed
    # in the previous iteration because it contained an incomplete
    # FASTA or FASTQ record.

    while True:
        if start == len(buf):
            raise OverflowError('FASTA/FASTQ record does not fit into buffer')
        bufend = f.readinto(memoryview(buf)[start:]) + start
        if start == bufend:
            # End of file
            break
        end = head(buf, bufend)
        assert end <= bufend
        if end > 0:
            yield memoryview(buf)[0:end]
        start = bufend - end
        assert start >= 0
        buf[0:start] = buf[end:bufend]

    if start > 0:
        yield memoryview(buf)[0:start]


def read_paired_chunks(f, f2, buffer_size=4*1024**2):
    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)

    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])
    if (start1 == 1 and buf1[0:1] != b'@') or (start2 == 1 and buf2[0:1] != b'@'):
        raise FileFormatError('Paired-end data must be in FASTQ format when using multiple cores')

    while True:
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = _paired_fastq_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0:
            yield (memoryview(buf1)[0:end1], memoryview(buf2)[0:end2])
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]

    if start1 > 0 or start2 > 0:
        yield (memoryview(buf1)[0:start1], memoryview(buf2)[0:start2])
