"""
Sequence I/O: Read and write FASTA and FASTQ files efficiently
"""

__all__ = [
    'open',
    'Sequence',
    'FastaReader',
    'FastaWriter',
    'FastqReader',
    'FastqWriter',
    'UnknownFileFormat',
    'FileFormatError',
    'FastaFormatError',
    'FastqFormatError',
    'InterleavedSequenceReader',
    'InterleavedSequenceWriter',
    'PairedSequenceReader',
]

import os
from contextlib import ExitStack
import functools

from xopen import xopen

from ._core import Sequence
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .exceptions import UnknownFileFormat, FileFormatError, FastaFormatError, FastqFormatError
from .chunks import read_chunks, read_paired_chunks

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def open(file1, file2=None, fileformat=None, interleaved=False, mode='r', qualities=None):
    """
    Open sequence files in FASTA or FASTQ format for reading or writing. This is
    a factory that returns an instance of one of the ...Reader or ...Writer
    classes also defined in this module.

    file1, file2 -- Paths to regular or compressed files or file-like
        objects. Use file1 if data is single-end. If also file2 is provided,
        sequences are paired.

    mode -- Either 'r' for reading or 'w' for writing.

    interleaved -- If True, then file1 contains interleaved paired-end data.
        file2 must be None in this case.

    fileformat -- If set to None, the file format is autodetected from the file name
        extension. Set to 'fasta' or 'fastq' to not auto-detect.

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
    if interleaved and file2 is not None:
        raise ValueError("When interleaved is set, file2 must be None")
    if file2 is not None:
        if mode == 'r':
            return PairedSequenceReader(file1, file2, fileformat)
        else:
            return PairedSequenceWriter(file1, file2, fileformat, qualities)
    if interleaved:
        if mode == 'r':
            return InterleavedSequenceReader(file1, fileformat)
        else:
            return InterleavedSequenceWriter(file1, fileformat, qualities)

    # The multi-file options have been dealt with, delegate rest to the
    # single-file function.
    return _open_single(
        file1, fileformat=fileformat, mode=mode, qualities=qualities)


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
    name, ext = os.path.splitext(name)
    if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
        return 'fasta'
    elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
        return 'fastq'
    return None


def _open_single(file, fileformat=None, mode='r', qualities=None):
    """
    Open a single sequence file. See description of open() above.
    """
    if mode not in ('r', 'w'):
        raise ValueError("Mode must be 'r' or 'w'")
    if isinstance(file, str):
        file = xopen(file, mode + 'b')
        close_file = True
    else:
        if not hasattr(file, 'readinto'):
            raise ValueError(
                'When passing in an open file-like object, it must have been opened in binary mode')
        close_file = False
    if mode == 'r':
        fastq_handler = FastqReader
        fasta_handler = FastaReader
    else:
        fastq_handler = FastqWriter
        fasta_handler = FastaWriter
    fastq_handler = functools.partial(fastq_handler, _close_file=close_file)
    fasta_handler = functools.partial(fasta_handler, _close_file=close_file)

    if fileformat:  # Explict file format given
        fileformat = fileformat.lower()
        if fileformat == 'fasta':
            return fasta_handler(file)
        elif fileformat == 'fastq':
            return fastq_handler(file)
        else:
            raise UnknownFileFormat(
                "File format {!r} is unknown (expected 'fasta' or 'fastq').".format(fileformat))

    # First, try to detect the file format from the file name only
    format = None
    if hasattr(file, "name") and isinstance(file.name, str):
        format = _detect_format_from_name(file.name)
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
            raise UnknownFileFormat(
                'Could not determine whether file {!r} is FASTA or FASTQ. The file extension was '
                'not available or not recognized and the first character in the file ({!r}) is '
                'unexpected.'.format(file, first_char))

    if format is None:
        assert mode == 'w'
        raise UnknownFileFormat('Cannot determine whether to write in FASTA or FASTQ format')

    if format == 'fastq' and mode == 'w' and qualities is False:
        raise ValueError(
            'Output format cannot be FASTQ since no quality values are available.')

    return fastq_handler(file) if format == 'fastq' else fasta_handler(file)


def _sequence_names_match(r1, r2):
    """
    Check whether the sequence records r1 and r2 have identical names, ignoring a
    suffix of '1' or '2'. Some old paired-end reads have names that end in '/1'
    and '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ)
    appends a .1 and .2 to paired-end reads if option -I is used.
    """
    name1 = r1.name.split(None, 1)[0]
    name2 = r2.name.split(None, 1)[0]
    if name1[-1:] in '12' and name2[-1:] in '12':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2


class PairedSequenceReader:
    """
    Read paired-end reads from two files.

    Wraps two BinaryFileReader instances, making sure that reads are properly
    paired.
    """
    paired = True

    def __init__(self, file1, file2, fileformat=None):
        with ExitStack() as stack:
            self.reader1 = stack.enter_context(_open_single(file1, fileformat=fileformat))
            self.reader2 = stack.enter_context(_open_single(file2, fileformat=fileformat))
            self._close = stack.pop_all().close
        self.delivers_qualities = self.reader1.delivers_qualities

    def __iter__(self):
        """
        Iterate over the paired reads. Each item is a pair of Sequence objects.
        """
        # Avoid usage of zip() below since it will consume one item too many.
        it1, it2 = iter(self.reader1), iter(self.reader2)
        while True:
            try:
                r1 = next(it1)
            except StopIteration:
                # End of file 1. Make sure that file 2 is also at end.
                try:
                    next(it2)
                    raise FileFormatError("Reads are improperly paired. There are more reads in "
                        "file 2 than in file 1.", line=None) from None
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FileFormatError("Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.", line=None) from None
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Read name '{}' "
                    "in file 1 does not match '{}' in file 2.".format(r1.name, r2.name), line=None) from None
            yield (r1, r2)

    def close(self):
        self._close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class InterleavedSequenceReader:
    """
    Read paired-end reads from an interleaved FASTQ file.
    """
    paired = True

    def __init__(self, file, fileformat=None):
        self.reader = _open_single(file, fileformat=fileformat)
        self.delivers_qualities = self.reader.delivers_qualities

    def __iter__(self):
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FileFormatError("Interleaved input file incomplete: Last record "
                    "{!r} has no partner.".format(r1.name), line=None) from None
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Name {!r} "
                    "(first) does not match {!r} (second).".format(r1.name, r2.name), line=None)
            yield (r1, r2)

    def close(self):
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class PairedSequenceWriter:
    def __init__(self, file1, file2, fileformat='fastq', qualities=None):
        with ExitStack() as stack:
            self._writer1 = stack.enter_context(_open_single(file1, fileformat=fileformat, mode='w',
                qualities=qualities))
            self._writer2 = stack.enter_context(_open_single(file2, fileformat=fileformat, mode='w',
                qualities=qualities))
            self._close = stack.pop_all().close

    def write(self, read1, read2):
        self._writer1.write(read1)
        self._writer2.write(read2)

    def close(self):
        self._close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceWriter:
    """
    Write paired-end reads to an interleaved FASTA or FASTQ file
    """

    def __init__(self, file, fileformat='fastq', qualities=None):

        self._writer = _open_single(
            file, fileformat=fileformat, mode='w', qualities=qualities)

    def write(self, read1, read2):
        self._writer.write(read1)
        self._writer.write(read2)

    def close(self):
        self._writer.close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()
