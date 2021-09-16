import functools
import os
from typing import Optional, Union, BinaryIO

from .exceptions import UnknownFileFormat
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from ._util import _is_path


def _open_single(
    file_or_path: Union[str, os.PathLike, BinaryIO],
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
        path = os.fspath(file_or_path)  # type: ignore
        file = opener(path, mode + "b")
        close_file = True
    else:
        if mode == "r" and not hasattr(file_or_path, "readinto"):
            raise ValueError(
                "When passing in an open file-like object, it must have been opened in binary mode"
            )
        file = file_or_path
        if hasattr(file, "name") and isinstance(file.name, str):
            path = file.name
        else:
            path = None
        close_file = False
    fastq_handler = FastqReader if mode == "r" else FastqWriter
    fasta_handler = FastaReader if mode == "r" else FastaWriter
    handlers = {
        "fastq": functools.partial(fastq_handler, _close_file=close_file),
        "fasta": functools.partial(fasta_handler, _close_file=close_file),
    }
    if fileformat:
        try:
            handler = handlers[fileformat.lower()]
        except KeyError:
            file.close()
            raise UnknownFileFormat(
                f"File format '{fileformat}' is unknown (expected 'fasta' or 'fastq')."
            )
        return handler(file)

    if path is not None:
        fileformat = _detect_format_from_name(path)
    if fileformat is None and mode == "w" and qualities is not None:
        # Format not recognized, but we know whether to use a format with or without qualities
        fileformat = "fastq" if qualities else "fasta"

    if mode == "r" and fileformat is None:
        fileformat = _detect_format_from_content(file)
        if fileformat is None:
            name = getattr(file, "name", repr(file))
            file.close()
            raise UnknownFileFormat(
                f"Could not determine whether file '{name}' is FASTA or FASTQ. The file extension "
                "is not available or not recognized, and the first character in the file is "
                "unexpected."
            )

    if fileformat is None:
        assert mode == "w"
        extra = " because the output file name is not available" if path is None else ""
        file.close()
        raise UnknownFileFormat(
            "Auto-detection of the output file format (FASTA/FASTQ) failed" + extra
        )

    if fileformat == "fastq" and mode in "wa" and qualities is False:
        file.close()
        raise ValueError(
            "Output format cannot be FASTQ since no quality values are available."
        )

    return handlers[fileformat](file)


def _detect_format_from_name(name: str) -> Optional[str]:
    """
    name -- file name

    Return 'fasta', 'fastq' or None if the format could not be detected.
    """
    name = name.lower()
    for ext in (".gz", ".xz", ".bz2"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    name, ext = os.path.splitext(name)
    if ext in [".fasta", ".fa", ".fna", ".csfasta", ".csfa"]:
        return "fasta"
    elif ext in [".fastq", ".fq"] or (ext == ".txt" and name.endswith("_sequence")):
        return "fastq"
    return None


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
        b"@": "fastq",
        b">": "fasta",
        b"#": "fasta",  # Some FASTA variants allow comments
        b"": "fastq",  # Pretend FASTQ for empty input
    }
    return formats.get(first_char, None)
