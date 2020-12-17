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
    'read_chunks',
    'read_paired_chunks',
    'record_names_match',
    '__version__',
]

import os
from os import fspath, PathLike
from contextlib import ExitStack
import functools
from typing import Optional, Union, BinaryIO, Tuple, Iterator

from xopen import xopen

from ._core import Sequence, record_names_match
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .exceptions import UnknownFileFormat, FileFormatError, FastaFormatError, FastqFormatError
from .chunks import read_chunks, read_paired_chunks
from ._version import version as __version__
from ._util import _is_path


def open(
    file1: Union[str, PathLike, BinaryIO],
    *,
    file2: Optional[Union[str, PathLike, BinaryIO]] = None,
    fileformat: Optional[str] = None,
    interleaved: bool = False,
    mode: str = "r",
    qualities: Optional[bool] = None,
    opener=xopen
) -> Union[
    FastaReader,
    FastaWriter,
    FastqReader,
    FastqWriter,
    "PairedSequenceReader",
    "PairedSequenceWriter",
    "PairedSequenceAppender",
    "InterleavedSequenceReader",
    "InterleavedSequenceWriter",
    "InterleavedSequenceAppender",
]:
    """
    Open sequence files in FASTA or FASTQ format for reading or writing. This is
    a factory that returns an instance of one of the ...Reader or ...Writer
    classes also defined in this module.

    file1, file2 -- Paths to regular or compressed files or file-like
        objects (as str or as pathlib.Path). Use only file1 if data is single-end.
        If sequences are paired, use also file2.

    mode -- Either 'r' for reading, 'w' for writing or 'a' for appending.

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

    opener -- A function that is used to open file1 and file2 if they are not
        already open file-like objects. By default, xopen is used, which can
        also open compressed file formats.
    """
    if mode not in ("r", "w", "a"):
        raise ValueError("Mode must be 'r', 'w' or 'a'")
    if interleaved and file2 is not None:
        raise ValueError("When interleaved is set, file2 must be None")

    if file2 is not None:
        if mode in "wa" and file1 == file2:
            raise ValueError("The paired-end output files are identical")
        if mode == "r":
            return PairedSequenceReader(file1, file2, fileformat, opener=opener)
        elif mode == "w":
            return PairedSequenceWriter(file1, file2, fileformat, qualities, opener=opener)
        else:
            return PairedSequenceAppender(file1, file2, fileformat, qualities, opener=opener)
    if interleaved:
        if mode == "r":
            return InterleavedSequenceReader(file1, fileformat, opener=opener)
        elif mode == "w":
            return InterleavedSequenceWriter(file1, fileformat, qualities, opener=opener)
        else:
            return InterleavedSequenceAppender(file1, fileformat, qualities, opener=opener)

    # The multi-file options have been dealt with, delegate rest to the
    # single-file function.
    return _open_single(
        file1, opener=opener, fileformat=fileformat, mode=mode, qualities=qualities)


def _detect_format_from_name(name: str) -> Optional[str]:
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


def _open_single(
    file_or_path: Union[str, PathLike, BinaryIO],
    opener,
    *,
    fileformat: Optional[str] = None,
    mode: str = "r",
    qualities: Optional[bool] = None,
) -> Union[FastaReader, FastaWriter, FastqReader, FastqWriter]:
    """
    Open a single sequence file. See description of open() above.
    """
    if mode not in ("r", "w", "a"):
        raise ValueError("Mode must be 'r', 'w' or 'a'")

    path: Optional[str]
    if _is_path(file_or_path):
        path = fspath(file_or_path)  # type: ignore
        file = opener(path, mode + "b")
        close_file = True
    else:
        if mode == 'r' and not hasattr(file_or_path, 'readinto'):
            raise ValueError(
                'When passing in an open file-like object, it must have been opened in binary mode')
        file = file_or_path
        if hasattr(file, "name") and isinstance(file.name, str):
            path = file.name
        else:
            path = None
        close_file = False
    fastq_handler = FastqReader if mode == "r" else FastqWriter
    fasta_handler = FastaReader if mode == "r" else FastaWriter
    handlers = {
        'fastq': functools.partial(fastq_handler, _close_file=close_file),
        'fasta': functools.partial(fasta_handler, _close_file=close_file),
    }
    if fileformat:
        try:
            handler = handlers[fileformat.lower()]
        except KeyError:
            file.close()
            raise UnknownFileFormat(
                "File format {!r} is unknown (expected 'fasta' or 'fastq').".format(fileformat))
        return handler(file)

    if path is not None:
        fileformat = _detect_format_from_name(path)
    if fileformat is None and mode == 'w' and qualities is not None:
        # Format not recognized, but we know whether to use a format with or without qualities
        fileformat = 'fastq' if qualities else 'fasta'

    if mode == 'r' and fileformat is None:
        fileformat = _detect_format_from_content(file)
        if fileformat is None:
            name = getattr(file, "name", repr(file))
            file.close()
            raise UnknownFileFormat(
                'Could not determine whether file "{}" is FASTA or FASTQ. The file extension was '
                'not available or not recognized and the first character in the file is '
                'unexpected.'.format(name))

    if fileformat is None:
        assert mode == 'w'
        extra = " because the output file name is not available" if path is None else ""
        file.close()
        raise UnknownFileFormat(
            "Auto-detection of the output file format (FASTA/FASTQ) failed" + extra)

    if fileformat == 'fastq' and mode in "wa" and qualities is False:
        file.close()
        raise ValueError(
            'Output format cannot be FASTQ since no quality values are available.')

    return handlers[fileformat](file)


def _detect_format_from_content(file: BinaryIO) -> Optional[str]:
    """
    Return 'fasta', 'fastq' or None
    """
    if file.seekable():
        first_char = file.read(1)
        if file.tell() > 0:
            file.seek(-1, 1)
    else:
        first_char = file.peek(1)[0:1]  # type: ignore
    formats = {
        b'@': 'fastq',
        b'>': 'fasta',
        b'#': 'fasta',  # Some FASTA variants allow comments
        b'': 'fastq',  # Pretend FASTQ for empty input
    }
    return formats.get(first_char, None)


class PairedSequenceReader:
    """
    Read paired-end reads from two files.

    Wraps two BinaryFileReader instances, making sure that reads are properly
    paired.
    """
    paired = True

    def __init__(
        self,
        file1: Union[str, PathLike, BinaryIO],
        file2: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        with ExitStack() as stack:
            self.reader1 = stack.enter_context(_open_single(file1, opener=opener, fileformat=fileformat))
            self.reader2 = stack.enter_context(_open_single(file2, opener=opener, fileformat=fileformat))
            self._close = stack.pop_all().close
        self.delivers_qualities = self.reader1.delivers_qualities

    def __repr__(self) -> str:
        return "PairedSequenceReader(file1={}, file2={})".format(self.reader1, self.reader2)

    def __iter__(self) -> Iterator[Tuple[Sequence, Sequence]]:
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
                    raise FileFormatError(
                        "Reads are improperly paired. There are more reads in "
                        "file 2 than in file 1.", line=None) from None
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FileFormatError(
                    "Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.", line=None) from None
            if not record_names_match(r1.name, r2.name):
                raise FileFormatError(
                    "Reads are improperly paired. Read name '{}' "
                    "in file 1 does not match '{}' in file 2.".format(r1.name, r2.name), line=None) from None
            yield (r1, r2)

    def close(self) -> None:
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

    def __init__(
        self,
        file: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        reader = _open_single(file, opener=opener, fileformat=fileformat)
        assert isinstance(reader, (FastaReader, FastqReader))  # for Mypy
        self.reader = reader
        self.delivers_qualities = self.reader.delivers_qualities

    def __repr__(self) -> str:
        return "InterleavedSequenceReader({})".format(self.reader)

    def __iter__(self) -> Iterator[Tuple[Sequence, Sequence]]:
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FileFormatError(
                    "Interleaved input file incomplete: Last record "
                    "{!r} has no partner.".format(r1.name), line=None) from None
            if not record_names_match(r1.name, r2.name):
                raise FileFormatError(
                    "Reads are improperly paired. Name {!r} "
                    "(first) does not match {!r} (second).".format(r1.name, r2.name), line=None)
            yield (r1, r2)

    def close(self) -> None:
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class PairedSequenceWriter:
    _mode = "w"

    def __init__(
        self,
        file1: Union[str, PathLike, BinaryIO],
        file2: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = "fastq",
        qualities: Optional[bool] = None,
        opener=xopen,
    ):
        with ExitStack() as stack:
            self._writer1: Union[FastaWriter, FastqWriter]
            self._writer2: Union[FastaWriter, FastqWriter]
            self._writer1 = stack.enter_context(
                _open_single(
                    file1, opener=opener, fileformat=fileformat, mode=self._mode, qualities=qualities))
            self._writer2 = stack.enter_context(
                _open_single(
                    file2, opener=opener, fileformat=fileformat, mode=self._mode, qualities=qualities))
            self._close = stack.pop_all().close

    def __repr__(self) -> str:
        return "{}({}, {})".format(self.__class__.__name__, self._writer1, self._writer2)

    def write(self, read1, read2) -> None:
        self._writer1.write(read1)
        self._writer2.write(read2)

    def close(self) -> None:
        self._close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()


class PairedSequenceAppender(PairedSequenceWriter):
    _mode = "a"


class InterleavedSequenceWriter:
    """
    Write paired-end reads to an interleaved FASTA or FASTQ file
    """
    _mode = "w"

    def __init__(
        self,
        file: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = "fastq",
        qualities: Optional[bool] = None,
        opener=xopen,
    ):
        writer = _open_single(
            file, opener=opener, fileformat=fileformat, mode=self._mode, qualities=qualities)
        assert isinstance(writer, (FastaWriter, FastqWriter))  # only for Mypy
        self._writer = writer

    def __repr__(self) -> str:
        return "{}({})".format(self.__class__.__name__, self._writer)

    def write(self, read1: Sequence, read2: Sequence) -> None:
        self._writer.write(read1)
        self._writer.write(read2)

    def close(self) -> None:
        self._writer.close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceAppender(InterleavedSequenceWriter):
    _mode = "a"
