"""
Chunked reading of FASTA and FASTQ files

This can be used to very quickly split up the input file into similarly-sized chunks,
without actually parsing the records. The chunks can then be distributed to worker threads
or subprocess and be parsed and processed there.
"""

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
    pos = buf.rfind(b"\n>", 0, end)
    if pos != -1:
        return pos + 1
    if buf[0:1] == b">" or buf[0:1] == b"#":
        return 0
    if len(buf) == 0:
        return 0
    c = chr(buf[0])
    raise FastaFormatError(
        f"FASTA file expected to start with '>', but found {repr(c)}",
        line=None,
    )


def _fastq_head(buf: bytes, end: Optional[int] = None) -> int:
    """
    Search for the end of the last complete *two* FASTQ records in buf[:end].

    Two FASTQ records are required to ensure that read pairs in interleaved
    paired-end data are not split.
    """
    linebreaks = buf.count(b"\n", 0, end)
    right = end
    for _ in range(linebreaks % 8 + 1):
        right = buf.rfind(b"\n", 0, right)
    # Note that this works even if linebreaks == 0:
    # rfind() returns -1 and adding 1 gives index 0,
    # which is correct.
    return right + 1  # type: ignore


def read_chunks(f: RawIOBase, buffer_size: int = 4 * 1024**2) -> Iterator[memoryview]:
    """
    Read chunks of complete FASTA or FASTQ records from a file.
    If the format is detected to be FASTQ, all chunks except possibly the last contain
    an even number of records such that interleaved paired-end reads remain in sync.
    The yielded memoryview objects are only valid for one iteration because the internal
    buffer is re-used in the next iteration.

    Arguments:
        f: File with FASTA or FASTQ reads; must have been opened in binary mode
        buffer_size: Largest allowed chunk size

    Yields:
        memoryview representing the chunk. This becomes invalid on the next iteration.

    Raises:
         ValueError: A FASTQ record was encountered that is larger than *buffer_size*.
         UnknownFileFormat: The file format could not be detected
           (the first byte must be "@", ">" or "#")
    """
    # This buffer is re-used in each iteration.
    buf = bytearray(buffer_size)

    # Read one byte to determine file format.
    # If there is a comment char, we assume FASTA!
    start = f.readinto(memoryview(buf)[0:1])
    if start == 0:
        # Empty file
        return
    assert start == 1
    if buf[0:1] == b"@":
        head = _fastq_head
    elif buf[0:1] == b"#" or buf[0:1] == b">":
        head = _fasta_head
    else:
        raise UnknownFileFormat(
            f"Cannnot determine input file format: First character expected to be '>' or '@', "
            f"but found {repr(chr(buf[0]))}"
        )

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
            raise OverflowError("FASTA/FASTQ record does not fit into buffer")
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
    """
    Read chunks of paired-end FASTQ reads from two files.
    A pair of chunks (memoryview objects) is yielded on each iteration,
    and both chunks are guaranteed to have the same number of sequences.
    That is, the paired-end reads will stay in sync.

    The memoryviews are only valid for one iteration because the internal
    buffer is re-used in the next iteration.

    This is similar to `read_chunks`, but for paired-end data.
    Unlike `read_chunks`, this only works for FASTQ input.

    Args:
        f: File with R1 reads; must have been opened in binary mode
        f2: File with R2 reads; must have been opened in binary mode
        buffer_size: Largest allowed chunk size

    Yields:
        Pairs of memoryview objects.

    Raises:
         ValueError: A FASTQ record was encountered that is larger than *buffer_size*.
    """
    if buffer_size < 6:
        raise ValueError("Buffer size too small")

    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)

    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])
    if (start1 == 1 and buf1[0:1] != b"@") or (start2 == 1 and buf2[0:1] != b"@"):
        raise FileFormatError(
            "Paired-end data must be in FASTQ format when using multiple cores",
            line=None,
        )

    while True:
        if start1 == len(buf1) and start2 == len(buf2):
            raise ValueError(
                f"FASTQ records do not fit into buffer of size {buffer_size}"
            )
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1  # type: ignore
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2  # type: ignore
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = _paired_fastq_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0:
            yield (memoryview(buf1)[0:end1], memoryview(buf2)[0:end2])
        else:
            assert end1 == 0 and end2 == 0
            extra = ""
            if bufend1 == 0 or bufend2 == 0:
                i = 1 if bufend1 == 0 else 2
                extra = f". File {i} ended, but more data found in the other file"
            raise FileFormatError(
                f"Premature end of paired FASTQ input{extra}.", line=None
            )
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]

    if start1 > 0 or start2 > 0:
        yield (memoryview(buf1)[0:start1], memoryview(buf2)[0:start2])
