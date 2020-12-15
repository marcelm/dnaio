"""Chunked reading of FASTA and FASTQ files"""
from io import RawIOBase
from typing import Optional, Iterator, Tuple

from ._core import paired_fastq_heads as _paired_fastq_heads
from .exceptions import FileFormatError, FastaFormatError, UnknownFileFormat


def _fasta_head(buf: bytes, end: Optional[int] = None) -> int:
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
    raise FastaFormatError('File does not start with ">"', line=None)


def _fastq_head(buf: bytes, end: Optional[int] = None) -> int:
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
    return right + 1  # type: ignore


def read_chunks(f: RawIOBase, buffer_size: int = 4 * 1024**2) -> Iterator[memoryview]:
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
    elif start == 0:
        # Empty file
        return
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
        bufend = f.readinto(memoryview(buf)[start:]) + start  # type: ignore
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


def read_paired_chunks(
    f: RawIOBase,
    f2: RawIOBase,
    buffer_size: int = 4 * 1024**2,
) -> Iterator[Tuple[memoryview, memoryview]]:
    if buffer_size < 1:
        raise ValueError("Buffer size too small")

    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)

    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])  # type: ignore
    start2 = f2.readinto(memoryview(buf2)[0:1])  # type: ignore
    if (start1 == 1 and buf1[0:1] != b'@') or (start2 == 1 and buf2[0:1] != b'@'):
        raise FileFormatError(
            "Paired-end data must be in FASTQ format when using multiple cores", line=None)

    while True:
        if start1 == len(buf1) or start2 == len(buf2):
            raise ValueError("FASTQ record does not fit into buffer")
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1  # type: ignore
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2  # type: ignore
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
