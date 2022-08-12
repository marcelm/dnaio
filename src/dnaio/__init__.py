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
    "read_chunks",
    "read_paired_chunks",
    "records_are_mates",
    "__version__",
]

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
from .multipleend import MultipleFileReader, MultipleFileWriter, _open_multiple
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


# Backwards compatibility aliases
Sequence = SequenceRecord


def open(
    file1: Union[str, PathLike, BinaryIO],
    *files: Union[str, PathLike, BinaryIO],
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
    MultipleFileReader,
    MultipleFileWriter,
]:
    """
    Open one or two files in FASTA or FASTQ format for reading or writing.

    Parameters:

      file1:
        Path or an open file-like object. For reading single-end reads, this is
        the only required argument.

      file2:
        Path or an open file-like object. When reading paired-end reads from
        two files, set this to the second file.

      mode:
        Either ``'r'``, ``'w'`` for writing or ``'a'`` for appending.

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

    Return:
       A subclass of `SingleEndReader`, `PairedEndReader`, `SingleEndWriter` or
       `PairedEndWriter`.
    """
    if mode not in ("r", "w", "a"):
        raise ValueError("Mode must be 'r', 'w' or 'a'")
    elif files and file2 is not None:
        raise ValueError(
            "the file2 argument can not be used when multiple "
            "input files are specified."
        )
    elif interleaved and (file2 is not None or files):
        raise ValueError(
            "When interleaved is True only one file must be specified as input."
        )
    if len(files) == 1:
        file2 = files[0]
    if interleaved or file2 is not None:
        return _open_paired(
            file1,
            file2=file2,
            opener=opener,
            fileformat=fileformat,
            interleaved=interleaved,
            mode=mode,
            qualities=qualities,
        )
    elif len(files) > 1:  # 3 or more files
        return _open_multiple(
            file1,
            *files,
            fileformat=fileformat,
            mode=mode,
            qualities=qualities,
            opener=opener
        )

    else:
        return _open_single(
            file1, opener=opener, fileformat=fileformat, mode=mode, qualities=qualities
        )
