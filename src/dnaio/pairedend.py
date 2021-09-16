from contextlib import ExitStack
from os import PathLike
from typing import Union, BinaryIO, Optional, Iterator, Tuple

from xopen import xopen

from ._core import Sequence, record_names_match
from .exceptions import FileFormatError
from .interfaces import PairedEndReader, PairedEndWriter
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .singleend import _open_single


class TwoFilePairedEndReader(PairedEndReader):
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
        *,
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        with ExitStack() as stack:
            self.reader1 = stack.enter_context(
                _open_single(file1, opener=opener, fileformat=fileformat)
            )
            self.reader2 = stack.enter_context(
                _open_single(file2, opener=opener, fileformat=fileformat)
            )
            self._close = stack.pop_all().close
        self.delivers_qualities = self.reader1.delivers_qualities

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(file1={self.reader1}, file2={self.reader2})"

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
                        "file 2 than in file 1.",
                        line=None,
                    ) from None
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FileFormatError(
                    "Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.",
                    line=None,
                ) from None
            if not record_names_match(r1.name, r2.name):
                raise FileFormatError(
                    f"Reads are improperly paired. Read name '{r1.name}' "
                    f"in file 1 does not match '{r2.name}' in file 2.",
                    line=None,
                ) from None
            yield (r1, r2)

    def close(self) -> None:
        self._close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class InterleavedPairedEndReader(PairedEndReader):
    """
    Read paired-end reads from an interleaved FASTQ file.
    """

    paired = True

    def __init__(
        self,
        file: Union[str, PathLike, BinaryIO],
        *,
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        reader = _open_single(file, opener=opener, fileformat=fileformat)
        assert isinstance(reader, (FastaReader, FastqReader))  # for Mypy
        self.reader = reader
        self.delivers_qualities = self.reader.delivers_qualities

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reader})"

    def __iter__(self) -> Iterator[Tuple[Sequence, Sequence]]:
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FileFormatError(
                    "Interleaved input file incomplete: Last record "
                    f"'{r1.name}' has no partner.",
                    line=None,
                ) from None
            if not record_names_match(r1.name, r2.name):
                raise FileFormatError(
                    f"Reads are improperly paired. Name '{r1.name}' "
                    f"(first) does not match '{r2.name}' (second).",
                    line=None,
                )
            yield (r1, r2)

    def close(self) -> None:
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class TwoFilePairedEndWriter(PairedEndWriter):
    def __init__(
        self,
        file1: Union[str, PathLike, BinaryIO],
        file2: Union[str, PathLike, BinaryIO],
        *,
        fileformat: Optional[str] = "fastq",
        qualities: Optional[bool] = None,
        opener=xopen,
        append: bool = False,
    ):
        mode = "a" if append else "w"
        with ExitStack() as stack:
            self._writer1: Union[FastaWriter, FastqWriter]
            self._writer2: Union[FastaWriter, FastqWriter]
            self._writer1 = stack.enter_context(
                _open_single(
                    file1,
                    opener=opener,
                    fileformat=fileformat,
                    mode=mode,
                    qualities=qualities,
                )
            )
            self._writer2 = stack.enter_context(
                _open_single(
                    file2,
                    opener=opener,
                    fileformat=fileformat,
                    mode=mode,
                    qualities=qualities,
                )
            )
            self._close = stack.pop_all().close

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._writer1}, {self._writer2})"

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


class InterleavedPairedEndWriter(PairedEndWriter):
    """
    Write paired-end reads to an interleaved FASTA or FASTQ file
    """

    def __init__(
        self,
        file: Union[str, PathLike, BinaryIO],
        *,
        fileformat: Optional[str] = "fastq",
        qualities: Optional[bool] = None,
        opener=xopen,
        append: bool = False,
    ):
        mode = "a" if append else "w"
        writer = _open_single(
            file, opener=opener, fileformat=fileformat, mode=mode, qualities=qualities
        )
        assert isinstance(writer, (FastaWriter, FastqWriter))  # only for Mypy
        self._writer = writer

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._writer})"

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
