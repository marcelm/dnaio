"""
Classes for reading FASTA and FASTQ files
"""
__all__ = ['FastaReader', 'FastqReader']

import io
from typing import Union, BinaryIO, Optional, Iterator, List

from xopen import xopen
from ._core import fastq_iter as _fastq_iter, Sequence
from ._util import shorten as _shorten
from .exceptions import FastaFormatError


class BinaryFileReader:
    """
    A mixin for readers that ensures that a file or a path can be passed in to the constructor.
    """
    _close_on_exit = False
    paired: bool = False
    mode: str = 'rb'

    def __init__(self, file: Union[str, BinaryIO], opener=xopen, _close_file: Optional[bool] = None):
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
        return "{}({!r})".format(self.__class__.__name__, getattr(self._file, "name", self._file))

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


class FastaReader(BinaryFileReader):
    """
    Reader for FASTA files.
    """

    def __init__(
        self,
        file: Union[str, BinaryIO],
        keep_linebreaks: bool = False,
        sequence_class=Sequence,
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
        self._delimiter = '\n' if keep_linebreaks else ''

    def __iter__(self) -> Iterator[Sequence]:
        """
        Read next entry from the file (single entry at a time).
        """
        name = None
        seq: List[str] = []
        f = io.TextIOWrapper(self._file)
        for i, line in enumerate(f):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == '>':
                if name is not None:
                    yield self.sequence_class(name, self._delimiter.join(seq), None)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FastaFormatError(
                    "Expected '>' at beginning of record, but got {!r}."
                    .format(_shorten(line)), line=i)

        if name is not None:
            yield self.sequence_class(name, self._delimiter.join(seq), None)
        # Prevent TextIOWrapper from closing the underlying file
        f.detach()


class FastqReader(BinaryFileReader):
    """
    Reader for FASTQ files. Does not support multi-line FASTQ files.
    """

    def __init__(
        self,
        file: Union[str, BinaryIO],
        sequence_class=Sequence,
        buffer_size: int = 1048576,
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
        # The first value yielded by _fastq_iter indicates
        # whether the file has repeated headers
        self._iter = _fastq_iter(self._file, self.sequence_class, self.buffer_size)
        try:
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

    def __iter__(self) -> Iterator[Sequence]:
        return self._iter
