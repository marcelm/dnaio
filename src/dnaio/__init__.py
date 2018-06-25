"""
Sequence I/O: Read and write FASTA and FASTQ files efficiently
"""
__version__ = '0.1'

import sys
from os.path import splitext

from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .exceptions import UnknownFileType, FileFormatError
from .colorspace import (ColorspaceFastaReader, ColorspaceFastaQualReader, FastaQualReader,
    ColorspaceFastaQualReader, SRAColorspaceFastqReader, ColorspaceFastqReader,
    ColorspaceFastqWriter, ColorspaceFastaWriter)


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
    return _open_single(
        file1, colorspace=colorspace, fileformat=fileformat, mode=mode, qualities=qualities)


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
            raise UnknownFileType(
                "File format {0!r} is unknown (expected "
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

    def __init__(self, file1, file2, colorspace=False, fileformat=None):
        self.reader1 = open(file1, colorspace=colorspace, fileformat=fileformat)
        self.reader2 = open(file2, colorspace=colorspace, fileformat=fileformat)
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
                        "file 2 than in file 1.")
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FileFormatError("Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.")
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Read name '{0}' "
                    "in file 1 does not match '{1}' in file 2.".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader1.close()
        self.reader2.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceReader:
    """
    Read paired-end reads from an interleaved FASTQ file.
    """
    paired = True

    def __init__(self, file, colorspace=False, fileformat=None):
        self.reader = open(file, colorspace=colorspace, fileformat=fileformat)
        self.delivers_qualities = self.reader.delivers_qualities

    def __iter__(self):
        # Avoid usage of zip() below since it will consume one item too many.
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FileFormatError("Interleaved input file incomplete: Last record "
                    "{!r} has no partner.".format(r1.name))
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Name {0!r} "
                    "(first) does not match {1!r} (second).".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class PairedSequenceWriter:
    def __init__(self, file1, file2, colorspace=False, fileformat='fastq', qualities=None):

        self._writer1 = open(file1, colorspace=colorspace, fileformat=fileformat, mode='w',
             qualities=qualities)
        self._writer2 = open(file2, colorspace=colorspace, fileformat=fileformat, mode='w',
             qualities=qualities)

    def write(self, read1, read2):
        self._writer1.write(read1)
        self._writer2.write(read2)

    def close(self):
        self._writer1.close()
        self._writer2.close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceWriter:
    """
    Write paired-end reads to an interleaved FASTA or FASTQ file
    """

    def __init__(self, file, colorspace=False, fileformat='fastq', qualities=None):
        from . import open as dnaio_open  # import locally to avoid circular import

        self._writer = dnaio_open(
            file, colorspace=colorspace, fileformat=fileformat, mode='w', qualities=qualities)

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
