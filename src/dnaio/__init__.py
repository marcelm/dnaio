"""
Sequence I/O: Read and write FASTA and FASTQ files efficiently
"""

__all__ = [
    "open",
    "SequenceRecord",
    "SingleEndReader",
    "PairedEndReader",
    "SingleEndWriter",
    "PairedEndWriter",
    "FastaReader",
    "FastaWriter",
    "FastqReader",
    "FastqWriter",
    "UnknownFileFormat",
    "FileFormatError",
    "FastaFormatError",
    "FastqFormatError",
    "InterleavedPairedEndReader",
    "InterleavedPairedEndWriter",
    "TwoFilePairedEndReader",
    "TwoFilePairedEndWriter",
    "MultipleFileReader",
    "MultipleFastaWriter",
    "MultipleFastqWriter",
    "read_chunks",
    "read_paired_chunks",
    "records_are_mates",
    "__version__",
]

import functools
from os import PathLike
from typing import Optional, Union, BinaryIO

from xopen import xopen

from ._core import (
    SequenceRecord,
)
from ._core import record_names_match  # noqa: F401  # deprecated
from ._core import records_are_mates
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .singleend import _open_single
from .pairedend import (
    _open_paired,
    TwoFilePairedEndReader,
    TwoFilePairedEndWriter,
    InterleavedPairedEndReader,
    InterleavedPairedEndWriter,
)
from .multipleend import (
    MultipleFastaWriter,
    MultipleFastqWriter,
    MultipleFileReader,
    MultipleFileWriter,
    _open_multiple,
)
from .exceptions import (
    UnknownFileFormat,
    FileFormatError,
    FastaFormatError,
    FastqFormatError,
)
from .interfaces import (
    SingleEndReader,
    PairedEndReader,
    SingleEndWriter,
    PairedEndWriter,
)
from .chunks import read_chunks, read_paired_chunks
from ._version import version as __version__


# Backwards compatibility alias
Sequence = SequenceRecord


def open(
    *files: Union[str, PathLike, BinaryIO],
    file1: Optional[Union[str, PathLike, BinaryIO]] = None,
    file2: Optional[Union[str, PathLike, BinaryIO]] = None,
    fileformat: Optional[str] = None,
    interleaved: bool = False,
    mode: str = "r",
    qualities: Optional[bool] = None,
    opener=xopen,
    compression_level: int = 1,
    open_threads: int = 0,
) -> Union[
    SingleEndReader,
    PairedEndReader,
    SingleEndWriter,
    PairedEndWriter,
    MultipleFileReader,
    MultipleFileWriter,
]:
    """
    Open one or two files in FASTA or FASTQ format for reading or writing.

    Parameters:
      files:
        one or more Path or open file-like objects. One for single-end
        reads, two for paired-end reads etc. More than two files are also
        supported. At least one file is required.

      file1:
        Deprecated keyword argument for the first file.

      file2:
        Deprecated keyword argument for the second file.

      mode:
        Set to ``'r'`` for reading, ``'w'`` for writing or ``'a'`` for appending.

      interleaved:
        If True, then file1 contains interleaved paired-end data.
        file2 must be None in this case.

      fileformat:
        If *None*, the file format is autodetected from the file name
        extension. Set to ``'fasta'`` or ``'fastq'`` to not auto-detect.

      qualities:
        When mode is ``'w'`` and fileformat is *None*, this can be set
        to *True* or *False* to specify whether the written sequences will have
        quality values. This is used in two ways:

        - If the output format cannot be determined (unrecognized extension
          etc.), no exception is raised, but FASTA or FASTQ format is chosen
          appropriately.

        - When False (no qualities available), an exception is raised when the
          auto-detected output format is FASTQ.

      opener: A function that is used to open file1 and file2 if they are not
        already open file-like objects. By default, ``xopen`` is used, which can
        also open compressed file formats.

      open_threads: By default dnaio opens files in the main thread,
        when threads is greater than 0 external processes are opened for
        compressing and decompressing files. This decreases wall clock time
        at the cost of a little extra overhead. This parameter does not work
        when a custom opener is set.

      compression_level: By default dnaio uses compression level 1 for writing
        gzipped files as this is the fastest. A higher level can be set using
        this parameter. This parameter does not work when a custom opener is
        set.

    Return:
       A subclass of `SingleEndReader`, `PairedEndReader`, `SingleEndWriter` or
       `PairedEndWriter`.
    """
    if files and (file1 is not None):
        raise ValueError(
            "The file1 keyword argument cannot be used together with files specified"
            "as positional arguments"
        )
    elif len(files) > 1 and file2 is not None:
        raise ValueError(
            "The file2 argument cannot be used together with more than one "
            "file specified as positional argument"
        )
    elif file1 is not None and file2 is not None and files:
        raise ValueError(
            "file1 and file2 arguments cannot be used together with files specified"
            "as positional arguments"
        )
    elif file1 is not None and file2 is not None:
        files = (file1, file2)
    elif file2 is not None and len(files) == 1:
        files = (files[0], file2)

    if len(files) > 1 and interleaved:
        raise ValueError("When interleaved is True, only one file must be specified.")
    elif mode not in ("r", "w", "a"):
        raise ValueError("Mode must be 'r', 'w' or 'a'")

    if opener == xopen:
        opener = functools.partial(
            xopen, threads=open_threads, compresslevel=compression_level
        )
    if interleaved or len(files) == 2:
        return _open_paired(
            *files,
            opener=opener,
            fileformat=fileformat,
            mode=mode,
            qualities=qualities,
        )
    elif len(files) > 2:
        return _open_multiple(
            *files, fileformat=fileformat, mode=mode, qualities=qualities, opener=opener
        )

    else:
        return _open_single(
            files[0],
            opener=opener,
            fileformat=fileformat,
            mode=mode,
            qualities=qualities,
        )
