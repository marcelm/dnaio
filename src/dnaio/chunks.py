"""
Chunked reading of FASTA and FASTQ files

This can be used to very quickly split up the input file into similarly-sized chunks,
without actually parsing the records. The chunks can then be distributed to worker threads
or subprocess and be parsed and processed there.
"""

from io import BufferedIOBase
from typing import Optional, Iterator, Tuple

from ._bam import read_bam_header_after_magic
from ._core import paired_fastq_heads as _paired_fastq_heads
from ._core import bam_head as _bam_head
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


def _paired_fasta_heads(
    buf1: bytes, buf2: bytes, end1: int, end2: int
) -> Tuple[int, int]:
    """
    Return positions pos1, pos2 where right1 <= end1 and right2 <= end2
    such that buf1[:pos1] and buf2[:pos2] contain the same number of complete FASTA
    records.
    """
    if end1 == 0 or end2 == 0:
        return (0, 0)

    if (end1 > 0 and buf1[:1] != b">") or (end2 > 0 and buf2[:1] != b">"):
        raise FastaFormatError("FASTA file expected to start with '>'", line=None)

    # Count complete records
    n_records1 = buf1.count(b"\n>", 0, end1)
    n_records2 = buf2.count(b"\n>", 0, end2)
    n_records = min(n_records1, n_records2)

    pos1 = pos2 = 0
    while n_records > 0:
        pos1 = buf1.find(b"\n>", pos1, end1) + 1
        pos2 = buf2.find(b"\n>", pos2, end2) + 1
        n_records -= 1

    return (pos1, pos2)


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


def read_chunks(
    f: BufferedIOBase, buffer_size: int = 4 * 1024**2
) -> Iterator[memoryview]:
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
    start = f.readinto(memoryview(buf)[0:4])
    if start == 0:
        # Empty file
        return
    assert start == 4
    if buf[0:1] == b"@":
        head = _fastq_head
    elif buf[0:1] == b"#" or buf[0:1] == b">":
        head = _fasta_head
    elif buf[0:4] == b"BAM\x01":
        head = _bam_head
        _ = read_bam_header_after_magic(f)
        start = 0  # Skip header and start at the records.
    else:
        raise UnknownFileFormat(
            f"Cannnot determine input file format: First characters expected "
            f"to be '>'. '@', or 'BAM\1', but found {repr(buf[0:4])}"
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


def read_paired_chunks(
    f: BufferedIOBase,
    f2: BufferedIOBase,
    buffer_size: int = 4 * 1024**2,
) -> Iterator[Tuple[memoryview, memoryview]]:
    """
    Read chunks of paired-end FASTA or FASTQ records from two files.
    A pair of chunks (memoryview objects) is yielded on each iteration,
    and both chunks are guaranteed to have the same number of sequences.
    That is, the paired-end reads will stay in sync.

    The memoryviews are only valid for one iteration because the internal
    buffer is re-used in the next iteration.

    This is similar to `read_chunks`, but for paired-end data.

    Args:
        f: File with R1 reads; must have been opened in binary mode
        f2: File with R2 reads; must have been opened in binary mode
        buffer_size: Largest allowed chunk size

    Yields:
        Pairs of memoryview objects.

    Raises:
         ValueError: A FASTA or FASTQ record was encountered that is larger than *buffer_size*.
    """
    if buffer_size < 6:
        raise ValueError("Buffer size too small")

    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)

    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])

    if start1 == 0 and start2 == 0:
        return

    if (start1 == 0) != (start2 == 0):
        i = 2 if start1 == 0 else 1
        raise FileFormatError(
            f"Paired-end reads not in sync: File with R{i} reads is empty and the other is not",
            line=None,
        )

    if buf1[:1] == b"@" != buf2[:1] == b"@":
        raise FileFormatError(
            "Paired-end data must be in FASTQ format when using multiple cores",
            line=None,
        )

    if buf1[:1] == b"@":
        file_format = "FASTQ"
        paired_heads = _paired_fastq_heads
    elif buf1[:1] == b">":
        file_format = "FASTA"
        paired_heads = _paired_fasta_heads
    else:
        raise FileFormatError(
            "First character in input file must be '@' (FASTQ) or '>' (FASTA), "
            f"but found {buf1[:1]}",
            line=None,
        )

    while True:
        if start1 == len(buf1) and start2 == len(buf2):
            raise ValueError(
                f"FASTA/FASTQ records do not fit into buffer of size {buffer_size}"
            )
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = paired_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0 or file_format == "FASTA":
            yield (memoryview(buf1)[0:end1], memoryview(buf2)[0:end2])
        else:
            assert end1 == 0 and end2 == 0
            extra = ""
            if bufend1 == 0 or bufend2 == 0:
                i = 1 if bufend1 == 0 else 2
                extra = f". File {i} ended, but more data found in the other file"
            raise FileFormatError(
                f"Premature end of paired-end input{extra}.", line=None
            )
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]

    if start1 > 0 or start2 > 0:
        yield (memoryview(buf1)[0:start1], memoryview(buf2)[0:start2])
