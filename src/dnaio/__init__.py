"""
Sequence I/O: Read and write FASTA and FASTQ files efficiently
"""
__version__ = '0.1'

import sys
from os.path import splitext

from .exceptions import UnknownFileType, FileFormatError
from .readers import FastaReader, PairedSequenceReader, InterleavedSequenceReader
from .writers import FastaWriter, FastqWriter, PairedSequenceWriter, InterleavedSequenceWriter
from .colorspace import (ColorspaceFastaReader, ColorspaceFastaQualReader, FastaQualReader, ColorspaceFastaQualReader,
    SRAColorspaceFastqReader, ColorspaceFastqReader, ColorspaceFastqWriter, ColorspaceFastaWriter)

from ._sequence import Sequence
from ._core import FastqReader, two_fastq_heads


def open(file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
        interleaved=False, mode='r', qualities=None):
    """
    Open sequence files in FASTA or FASTQ format for reading or writing. This is
    a factory that returns an instance of one of the ...Reader or ...Writer
    classes also defined in this module.

    file1, file2, qualfile -- Paths to regular or compressed files or file-like
        objects. Use file1 if data is single-end. If also file2 is provided,
        sequences are paired. If qualfile is given, then file1 must be a FASTA
        file and sequences are single-end. One of file2 and qualfile must always
        be None (no paired-end data is supported when reading qualfiles).

    mode -- Either 'r' for reading or 'w' for writing.

    interleaved -- If True, then file1 contains interleaved paired-end data.
        file2 and qualfile must be None in this case.

    colorspace -- If True, instances of the Colorspace... classes
        are returned.

    fileformat -- If set to None, file format is autodetected from the file name
        extension. Set to 'fasta', 'fastq', or 'sra-fastq' to not auto-detect.
        Colorspace is not auto-detected and must always be requested explicitly.

    qualities -- When mode is 'w' and fileformat is None, this can be set to
        True or False to specify whether the written sequences will have quality
        values. This is is used in two ways:
        * If the output format cannot be determined (unrecognized extension
          etc), no exception is raised, but fasta or fastq format is chosen
          appropriately.
        * When False (no qualities available), an exception is raised when the
          auto-detected output format is FASTQ.
    """
    if mode not in ('r', 'w'):
        raise ValueError("Mode must be 'r' or 'w'")
    if interleaved and (file2 is not None or qualfile is not None):
        raise ValueError("When interleaved is set, file2 and qualfile must be None")
    if file2 is not None and qualfile is not None:
        raise ValueError("Setting both file2 and qualfile is not supported")
    if file2 is not None:
        if mode == 'r':
            return PairedSequenceReader(file1, file2, colorspace, fileformat)
        else:
            return PairedSequenceWriter(file1, file2, colorspace, fileformat, qualities)

    if interleaved:
        if mode == 'r':
            return InterleavedSequenceReader(file1, colorspace, fileformat)
        else:
            return InterleavedSequenceWriter(file1, colorspace, fileformat, qualities)

    if qualfile is not None:
        if mode == 'w':
            raise NotImplementedError('Writing to csfasta/qual not supported')
        if colorspace:
            # read from .(CS)FASTA/.QUAL
            return ColorspaceFastaQualReader(file1, qualfile)
        else:
            return FastaQualReader(file1, qualfile)

    # All the multi-file things have been dealt with, delegate rest to the
    # single-file function.
    return _open_single(file1, colorspace=colorspace, fileformat=fileformat,
        mode=mode, qualities=qualities)


def _detect_format_from_name(name):
    """
    name -- file name

    Return 'fasta', 'fastq' or None if the format could not be detected.
    """
    name = name.lower()
    for ext in ('.gz', '.xz', '.bz2'):
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    name, ext = splitext(name)
    if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
        return 'fasta'
    elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
        return 'fastq'
    return None


def _open_single(file, colorspace=False, fileformat=None, mode='r', qualities=None):
    """
    Open a single sequence file. See description above.
    """
    if mode == 'r':
        fastq_handler = ColorspaceFastqReader if colorspace else FastqReader
        fasta_handler = ColorspaceFastaReader if colorspace else FastaReader
    elif mode == 'w':
        fastq_handler = ColorspaceFastqWriter if colorspace else FastqWriter
        fasta_handler = ColorspaceFastaWriter if colorspace else FastaWriter
    else:
        raise ValueError("Mode must be 'r' or 'w'")

    if fileformat:  # Explict file format given
        fileformat = fileformat.lower()
        if fileformat == 'fasta':
            return fasta_handler(file)
        elif fileformat == 'fastq':
            return fastq_handler(file)
        elif fileformat == 'sra-fastq' and colorspace:
            if mode == 'w':
                raise NotImplementedError('Writing to sra-fastq not supported')
            return SRAColorspaceFastqReader(file)
        else:
            raise UnknownFileType("File format {0!r} is unknown (expected "
                "'sra-fastq' (only for colorspace), 'fasta' or 'fastq').".format(fileformat))

    # Detect file format
    name = None
    if file == "-":
        file = sys.stdin.buffer if mode == 'r' else sys.stdout.buffer
    elif isinstance(file, str):
        name = file
    elif hasattr(file, "name"):  # seems to be an open file-like object
        name = file.name
        if not hasattr(file, 'readinto'):
            raise ValueError(
                'When passing in an open file-like object, it must have been opened in binary mode')

    format = _detect_format_from_name(name) if name else None

    if format is None and mode == 'w' and qualities is not None:
        # Format not recognized, but we know whether to use a format with or without qualities
        format = 'fastq' if qualities else 'fasta'

    if mode == 'r' and format is None:
        # No format detected so far. Try to read from the file.
        if file.seekable():
            first_char = file.read(1)
            file.seek(-1, 1)
        else:
            first_char = file.peek(1)[0:1]
        if first_char == b'#':
            # A comment char - only valid for some FASTA variants (csfasta)
            format = 'fasta'
        elif first_char == b'>':
            format = 'fasta'
        elif first_char == b'@':
            format = 'fastq'
        elif first_char == b'':
            # Empty input. Pretend this is FASTQ
            format = 'fastq'
        else:
            raise UnknownFileType(
                'Could not determine whether file {!r} is FASTA or FASTQ. The file extension was '
                'not available or not recognized and the first character in the file ({!r}) is '
                'unexpected.'.format(file, first_char))

    if format is None:
        assert mode == 'w'
        raise UnknownFileType('Cannot determine whether to write in FASTA or FASTQ format')

    if format == 'fastq' and mode == 'w' and qualities is False:
        raise ValueError(
            'Output format cannot be FASTQ since no quality values are available.')

    return fastq_handler(file) if format == 'fastq' else fasta_handler(file)


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
    elif start == 1 and buf[0:1] == b'#' or buf[0:1] == b'>':
        head = _fasta_head
    elif start > 0:
        raise UnknownFileType('Input file format unknown')

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

        end1, end2 = two_fastq_heads(buf1, buf2, bufend1, bufend2)
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
