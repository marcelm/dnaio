from os import PathLike
from typing import Union, BinaryIO, Optional

from xopen import xopen

from . import Sequence
from ._util import _is_path
from .interfaces import SingleEndWriter


class FileWriter:
    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        if _is_path(file):
            self._file = opener(file, "wb")
            self._close_on_exit = True
        else:
            self._file = file
            self._close_on_exit = bool(_close_file)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{getattr(self._file, 'name', self._file)}')"

    def close(self) -> None:
        if self._close_on_exit:
            self._file.close()

    def __enter__(self):
        if self._file.closed:
            raise ValueError("I/O operation on closed file")
        return self

    def __exit__(self, *args):
        self.close()


class FastaWriter(FileWriter, SingleEndWriter):
    """
    Write FASTA-formatted sequences to a file.
    """

    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        line_length: Optional[int] = None,
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        """
        If line_length is not None, the lines will
        be wrapped after line_length characters.
        """
        super().__init__(file, opener=opener, _close_file=_close_file)
        self.line_length = line_length if line_length != 0 else None

    def __repr__(self) -> str:
        return f"FastaWriter('{getattr(self._file, 'name', self._file)}')"

    def write(self, name_or_record, sequence: Optional[str] = None):
        """Write an entry to the the FASTA file.

        If only one parameter (name_or_record) is given, it must have
        attributes .name and .sequence, which are then used.
        Otherwise, the first parameter must be the name and the second
        the sequence.

        The effect is that you can write this:
        writer.write("name", "ACCAT")
        or
        writer.write(Sequence("name", "ACCAT"))
        """
        if sequence is None:
            name = name_or_record.name
            sequence = name_or_record.sequence
        else:
            name = name_or_record

        if self.line_length is not None:
            self._file.write(('>' + name + '\n').encode('ascii'))
            s = []
            for i in range(0, len(sequence), self.line_length):
                s.append(sequence[i:i + self.line_length] + '\n')
            self._file.write(''.join(s).encode('ascii'))
        else:
            text = '>' + name + '\n' + sequence + '\n'
            self._file.write(text.encode('ascii'))


class FastqWriter(FileWriter, SingleEndWriter):
    """
    Write sequences with qualities in FASTQ format.

    FASTQ files are formatted like this::

        @read name
        AACCGGTT
        +
        FF,:F,,F
    """
    file_mode = 'wb'

    def __init__(
        self,
        file: Union[PathLike, str, BinaryIO],
        *,
        two_headers: bool = False,
        opener=xopen,
        _close_file: Optional[bool] = None,
    ):
        super().__init__(file, opener=opener, _close_file=_close_file)
        self._two_headers = two_headers
        # setattr avoids a complaint from Mypy
        setattr(self, "write", self._write_two_headers if self._two_headers else self._write)

    def __repr__(self) -> str:
        return f"FastqWriter('{getattr(self._file, 'name', self._file)}')"

    def write(self, record: Sequence) -> None:
        """
        Dummy method to make it possible to instantiate this class.
        The correct write method is assigned in the constructor.
        """
        assert False

    def _write(self, record: Sequence) -> None:
        """
        Write a Sequence record to the FASTQ file.

        """
        self._file.write(record.fastq_bytes())

    def _write_two_headers(self, record: Sequence) -> None:
        """
        Write a Sequence record to the FASTQ file, repeating the header
        in the third line after the "+" .
        """
        self._file.write(record.fastq_bytes_two_headers())

    def writeseq(self, name: str, sequence: str, qualities: str) -> None:
        self._file.write(f"@{name:s}\n{sequence:s}\n+\n{qualities:s}\n".encode("ascii"))
