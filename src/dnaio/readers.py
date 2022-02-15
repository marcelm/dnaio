"""
Classes for reading FASTA and FASTQ files
"""
__all__ = ['FastaReader', 'FastqReader']

import io
from typing import Union, BinaryIO, Optional, Iterator, List

from xopen import xopen

from ._core import FastqIter, Sequence
from ._util import shorten as _shorten
from .exceptions import FastaFormatError
from .interfaces import SingleEndReader


class BinaryFileReader:
    """
    A mixin for readers that ensures that a file or a path can be passed in to the constructor.
    """
    _close_on_exit = False
    paired: bool = False
    mode: str = 'rb'

    def __init__(
        self, file: Union[str, BinaryIO], *, opener=xopen, _close_file: Optional[bool] = None
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
    Reader for FASTA files.
    """

    def __init__(
        self,
        file: Union[str, BinaryIO],
        *,
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
        self.number_of_records = 0

    def __iter__(self) -> Iterator[Sequence]:
        """
        Read next entry from the file (single entry at a time).
        """
        name = None
        seq: List[str] = []
        if self._file.closed:
            return
        f = io.TextIOWrapper(self._file)
        for i, line in enumerate(f):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == '>':
                if name is not None:
                    self.number_of_records += 1
                    yield self.sequence_class(name, self._delimiter.join(seq), None)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FastaFormatError(
                    f"Expected '>' at beginning of record, but got '{_shorten(line)}'.", line=i)

        if name is not None:
            self.number_of_records += 1
            yield self.sequence_class(name, self._delimiter.join(seq), None)
        # Prevent TextIOWrapper from closing the underlying file
        f.detach()


class FastqReader(BinaryFileReader, SingleEndReader):
    """
    Reader for FASTQ files. Does not support multi-line FASTQ files.
    """

    def __init__(
        self,
        file: Union[str, BinaryIO],
        *,
        sequence_class=Sequence,
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
        self._iter: FastqIter = FastqIter(self._file, self.sequence_class, self.buffer_size)
        self._two_headers = None

    def __iter__(self) -> Iterator[Sequence]:
        return self._iter

    @property
    def two_headers(self):
        if self._two_headers is not None:
            return self._two_headers

        if hasattr(self._file, "peek"):
            peek = self._file.peek
        elif self._file.seekable():
            def peek(n):
                original_pos = self._file.tell()
                chunk = self._file.read(n)
                self._file.seek(original_pos)
                return chunk
        else:
            raise io.UnsupportedOperation

        peek_size = 512
        while peek_size <= 2 * 1024 * 1024:
            chunk = peek(peek_size)
            lines = chunk.split(b"\n")
            for i, line in enumerate(lines):
                if i + 2 >= len(lines):
                    break
                # Condition below is only true for the start of a record. Not
                # for a quality line that starts with @.
                if line.startswith(b"@") and lines[i+2].startswith(b"+"):
                    self._two_headers = (lines[i+2] != b"+" and
                                         lines[i+2] != b"+\r")
                    return self._two_headers
            peek_size *= 2
        # Fallback to False. This file is likely not a valid FASTQ file. This
        # will be discovered during iteration and a proper error will be thrown
        # there.
        return False

    @property
    def number_of_records(self):
        try:
            return self._iter.number_of_records
        except AttributeError:
            return 0
