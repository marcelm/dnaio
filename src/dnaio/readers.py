"""
Classes for reading FASTA and FASTQ files
"""
__all__ = ["FastaReader", "FastqReader"]

import io
from os import PathLike
from typing import Union, BinaryIO, Optional, Iterator, List

from xopen import xopen

from ._core import FastqIter, SequenceRecord
from ._util import shorten as _shorten
from .exceptions import FastaFormatError
from .interfaces import SingleEndReader


class BinaryFileReader:
    """
    A mixin for readers that ensures that a file or a path can be passed in to the constructor.
    """

    _close_on_exit = False
    paired: bool = False
    mode: str = "rb"

    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        """
        The file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).
        """
        if isinstance(file, str):
            self._file = opener(file, self.mode)
            self._close_on_exit = True
        elif _close_file:
            self._close_on_exit = True
            self._file = file
        else:
            self._file = file

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{getattr(self._file, 'name', self._file)}')"

    def close(self) -> None:
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        if self._file is None:
            raise ValueError("I/O operation on closed BinaryFileReader")
        return self

    def __exit__(self, *args):
        self.close()


class FastaReader(BinaryFileReader, SingleEndReader):
    """
    Reader for FASTA files

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
    """

    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        keep_linebreaks: bool = False,
        sequence_class=SequenceRecord,
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        """
        file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).

        keep_linebreaks -- whether to keep newline characters in the sequence
        """
        super().__init__(file, opener=opener, _close_file=_close_file)
        self.sequence_class = sequence_class
        self.delivers_qualities = False
        self._delimiter = "\n" if keep_linebreaks else ""
        self.number_of_records = 0
        self._file = io.TextIOWrapper(self._file)

    def __iter__(self) -> Iterator[SequenceRecord]:
        """
        Iterate over the records in this FASTA file.
        """
        name = None
        seq: List[str] = []
        if self._file.closed:
            return
        for i, line in enumerate(self._file):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == ">":
                if name is not None:
                    self.number_of_records += 1
                    try:
                        yield self.sequence_class(name, self._delimiter.join(seq), None)
                    except ValueError as e:
                        raise FastaFormatError(
                            str(e)
                            + " (line number refers to record after the problematic one)",
                            line=i,
                        )
                name = line[1:]
                seq = []
            elif line and line[0] == "#":
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FastaFormatError(
                    f"Expected '>' at beginning of record, but got '{_shorten(line)}'.",
                    line=i,
                )

        if name is not None:
            self.number_of_records += 1
            try:
                yield self.sequence_class(name, self._delimiter.join(seq), None)
            except ValueError as e:
                raise FastaFormatError(str(e), line=None)


class FastqReader(BinaryFileReader, SingleEndReader):
    """
    Reader for FASTQ files. Does not support multi-line FASTQ files.

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
    """

    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        sequence_class=SequenceRecord,
        buffer_size: int = 128 * 1024,  # Buffer size used by cat, pigz etc.
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        """
        super().__init__(file, opener=opener, _close_file=_close_file)
        self.sequence_class = sequence_class
        self.delivers_qualities = True
        self.buffer_size = buffer_size
        try:
            self._iter: Iterator[SequenceRecord] = FastqIter(
                self._file, self.sequence_class, self.buffer_size
            )
        except Exception:
            self.close()
            raise
        try:
            # The first value yielded by FastqIter indicates
            # whether the file has repeated headers
            th = next(self._iter)
            assert isinstance(th, bool)
            self.two_headers: bool = th
        except StopIteration:
            # Empty file
            self.two_headers = False
            self._iter = iter(())
        except Exception:
            self.close()
            raise

    def __iter__(self) -> Iterator[SequenceRecord]:
        """Iterate over the records in this FASTQ file."""
        return self._iter

    @property
    def number_of_records(self):
        try:
            return self._iter.number_of_records
        except AttributeError:
            return 0
