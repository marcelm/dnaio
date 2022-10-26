from contextlib import ExitStack
from os import PathLike
from typing import Union, BinaryIO, Optional, Iterator, Tuple

from xopen import xopen

from ._core import SequenceRecord
from .exceptions import FileFormatError
from .interfaces import PairedEndReader, PairedEndWriter
from .readers import FastaReader, FastqReader
from .writers import FastaWriter, FastqWriter
from .singleend import _open_single


def _open_paired(
    *files: Union[str, PathLike, BinaryIO],
    fileformat: Optional[str] = None,
    mode: str = "r",
    qualities: Optional[bool] = None,
    opener=xopen,
) -> Union[PairedEndReader, PairedEndWriter]:
    """
    Open paired-end reads
    """
    if len(files) == 2:
        if mode in "wa" and files[0] == files[1]:
            raise ValueError("The paired-end output files are identical")
        if "r" in mode:
            return TwoFilePairedEndReader(
                *files, fileformat=fileformat, opener=opener, mode=mode
            )
        append = mode == "a"
        return TwoFilePairedEndWriter(
            *files,
            fileformat=fileformat,
            qualities=qualities,
            opener=opener,
            append=append,
        )
    elif len(files) == 1:
        if "r" in mode:
            return InterleavedPairedEndReader(
                files[0], fileformat=fileformat, opener=opener, mode=mode
            )
        append = mode == "a"
        return InterleavedPairedEndWriter(
            files[0],
            fileformat=fileformat,
            qualities=qualities,
            opener=opener,
            append=append,
        )
    raise ValueError("_open_paired must be called with one or two files.")


class TwoFilePairedEndReader(PairedEndReader):
    """
    Read paired-end reads from two files (not interleaved)

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
    """

    paired = True

    def __init__(
        self,
        file1: Union[str, PathLike, BinaryIO],
        file2: Union[str, PathLike, BinaryIO],
        *,
        mode="r",
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        self.mode = mode
        with ExitStack() as stack:
            self.reader1 = stack.enter_context(
                _open_single(file1, opener=opener, fileformat=fileformat, mode=mode)
            )
            self.reader2 = stack.enter_context(
                _open_single(file2, opener=opener, fileformat=fileformat, mode=mode)
            )
            self._close = stack.pop_all().close
        self.delivers_qualities = self.reader1.delivers_qualities

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(file1={self.reader1}, file2={self.reader2})"

    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
        """
        Iterate over the paired reads.
        Each yielded item is a pair of `SequenceRecord` objects.

        Raises a `FileFormatError` if reads are improperly paired.
        """
        for r1, r2 in zip(self.reader1, self.reader2):
            if not r1.is_mate(r2):
                raise FileFormatError(
                    f"Reads are improperly paired. Read name '{r1.name}' "
                    f"in file 1 does not match '{r2.name}' in file 2.",
                    line=None,
                ) from None
            yield r1, r2

        # Force consumption of another read to test if iterators are out of sync.
        try:
            next(iter(self.reader1))
        except StopIteration:
            pass
        try:
            next(iter(self.reader2))
        except StopIteration:
            pass
        if self.reader1.number_of_records < self.reader2.number_of_records:
            raise FileFormatError(
                "Reads are improperly paired. There are more reads in "
                "file 2 than in file 1.",
                line=None,
            ) from None
        if self.reader1.number_of_records > self.reader2.number_of_records:
            raise FileFormatError(
                "Reads are improperly paired. There are more reads in "
                "file 1 than in file 2.",
                line=None,
            ) from None

    def close(self) -> None:
        self._close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class InterleavedPairedEndReader(PairedEndReader):
    """
    Read paired-end reads from an interleaved FASTQ file

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
    """

    paired = True

    def __init__(
        self,
        file: Union[str, PathLike, BinaryIO],
        *,
        mode="r",
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        self.mode = mode
        reader = _open_single(file, opener=opener, mode=mode, fileformat=fileformat)
        assert isinstance(reader, (FastaReader, FastqReader))  # for Mypy
        self.reader = reader
        self.delivers_qualities = self.reader.delivers_qualities

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reader})"

    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
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
            if not r1.is_mate(r2):
                raise FileFormatError(
                    f"Reads are improperly paired. Name '{r1.name}' "
                    f"(first) does not match '{r2.name}' (second).",
                    line=None,
                )
            yield r1, r2

    def close(self) -> None:
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class TwoFilePairedEndWriter(PairedEndWriter):
    """
    Write paired-end reads to two files (not interleaved)

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
    """

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

    While this class can be instantiated directly, the recommended way is to
    use `dnaio.open` with appropriate arguments.
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

    def write(self, read1: SequenceRecord, read2: SequenceRecord) -> None:
        self._writer.write(read1)
        self._writer.write(read2)

    def close(self) -> None:
        self._writer.close()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()
