"""
Sequence I/O: Read and write FASTA and FASTQ files efficiently
"""

__all__ = [
    'open',
    'Sequence',
    'SingleEndReader',
    'PairedEndReader',
    'SingleEndWriter',
    'PairedEndWriter',
    'FastaReader',
    'FastaWriter',
    'FastqReader',
    'FastqWriter',
    'UnknownFileFormat',
    'FileFormatError',
    'FastaFormatError',
    'FastqFormatError',
    'InterleavedPairedEndReader',
    'InterleavedPairedEndWriter',
    'TwoFilePairedEndReader',
    'TwoFilePairedEndWriter',
    'read_chunks',
    'read_paired_chunks',
    'record_names_match',
    '__version__',
]

from os import PathLike
from typing import Optional, Union, BinaryIO

from xopen import xopen

from ._core import Sequence, record_names_match
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .singleend import _open_single
from .pairedend import (
    TwoFilePairedEndReader,
    TwoFilePairedEndWriter,
    InterleavedPairedEndReader,
    InterleavedPairedEndWriter,
)
from .exceptions import UnknownFileFormat, FileFormatError, FastaFormatError, FastqFormatError
from .interfaces import SingleEndReader, PairedEndReader, SingleEndWriter, PairedEndWriter
from .chunks import read_chunks, read_paired_chunks
from ._version import version as __version__


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
    SingleEndReader,
    PairedEndReader,
    SingleEndWriter,
    PairedEndWriter,
]:
    """
    Open sequence files in FASTA or FASTQ format for reading or writing.

    Parameters:
      file1:
      file2:
        Paths to regular or compressed files or file-like
        objects (as str or as pathlib.Path). Use only file1 if data is single-end.
        If sequences are paired, use also file2.
      mode:
        Either ``'r'`` for reading, ``'w'`` for writing or ``'a'`` for appending.
      interleaved:
        If True, then file1 contains interleaved paired-end data.
        file2 must be None in this case.
      fileformat:
        If *None*, the file format is autodetected from the file name
        extension. Set to ``'fasta'`` or ``'fastq'`` to not auto-detect.
      qualities:
        When mode is ``'w'`` and fileformat is *None*, this can be set
        to *True* or *False* to specify whether the written sequences will have
        quality values. This is is used in two ways:

        - If the output format cannot be determined (unrecognized extension
          etc), no exception is raised, but fasta or fastq format is chosen
          appropriately.

        - When False (no qualities available), an exception is raised when the
          auto-detected output format is FASTQ.
      opener: A function that is used to open file1 and file2 if they are not
        already open file-like objects. By default, ``xopen`` is used, which can
        also open compressed file formats.

    Return:
       An instance of one of the ...Reader or ...Writer classes
    """
    if mode not in ("r", "w", "a"):
        raise ValueError("Mode must be 'r', 'w' or 'a'")
    if interleaved and file2 is not None:
        raise ValueError("When interleaved is set, file2 must be None")

    if file2 is not None:
        if mode in "wa" and file1 == file2:
            raise ValueError("The paired-end output files are identical")
        if mode == "r":
            return TwoFilePairedEndReader(file1, file2, fileformat=fileformat, opener=opener)
        append = mode == "a"
        return TwoFilePairedEndWriter(
            file1, file2, fileformat=fileformat, qualities=qualities, opener=opener, append=append
        )
    if interleaved:
        if mode == "r":
            return InterleavedPairedEndReader(file1, fileformat=fileformat, opener=opener)
        append = mode == "a"
        return InterleavedPairedEndWriter(
            file1, fileformat=fileformat, qualities=qualities, opener=opener, append=append)

    # The multi-file options have been dealt with, delegate rest to the
    # single-file function.
    return _open_single(
        file1, opener=opener, fileformat=fileformat, mode=mode, qualities=qualities)
