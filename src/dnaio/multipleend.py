from os import PathLike
from typing import BinaryIO, IO, Iterable, Iterator, List, Optional, Tuple, Union

from xopen import xopen

from ._core import SequenceRecord, records_are_mates
from .exceptions import FileFormatError
from .interfaces import MultipleFileWriter
from .readers import FastaReader, FastqReader
from .singleend import _open_single
from .writers import FastaWriter, FastqWriter


def _open_multiple(
    *files: Union[str, PathLike, BinaryIO],
    fileformat: Optional[str] = None,
    mode: str = "r",
    qualities: Optional[bool] = None,
    opener=xopen,
):
    if mode not in ("r", "w", "a"):
        raise ValueError("mode must be one of 'r', 'w', 'a'")
    if mode == "r":
        return MultipleFileReader(*files, fileformat=fileformat, opener=opener)
    append = True if mode == "a" else False
    if fileformat == "fastq" or qualities or (fileformat is None and qualities is None):
        return MultipleFastqWriter(*files, opener=opener, append=append)
    return MultipleFastaWriter(*files, opener=opener, append=append)


class MultipleFileReader:
    def __init__(
        self,
        *files: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = None,
        opener=xopen,
    ):
        if len(files) < 1:
            raise ValueError("At least one file is required")
        self.files = files
        self.readers: List[Union[FastaReader, FastqReader]] = [
            _open_single(  # type: ignore
                file, opener=opener, fileformat=fileformat, mode="r"
            )
            for file in self.files
        ]
        self.delivers_qualities: bool = self.readers[0].delivers_qualities

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}"
            f"({', '.join(repr(reader) for reader in self.readers)})"
        )

    def __iter__(self) -> Iterator[Tuple[SequenceRecord, ...]]:
        if len(self.files) == 1:
            yield from zip(self.readers[0])
        else:
            for records in zip(*self.readers):
                if not records_are_mates(*records):
                    raise FileFormatError(
                        f"Records are out of sync, names "
                        f"{', '.join(repr(r.name) for r in records)} do not match.",
                        line=None,
                    )
                yield records
        # Consume one iteration to check if all the files have an equal number
        # of records..
        for reader in self.readers:
            try:
                _ = next(iter(reader))
            except StopIteration:
                pass
        record_numbers = [r.number_of_records for r in self.readers]
        if len(set(record_numbers)) != 1:
            raise FileFormatError(
                f"Files: {', '.join(str(file) for file in self.files)} have "
                f"an unequal amount of reads.",
                line=None,
            )

    def close(self):
        for reader in self.readers:
            reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class MultipleFastaWriter(MultipleFileWriter):
    def __init__(
        self,
        *files: Union[str, PathLike, BinaryIO],
        opener=xopen,
        append: bool = False,
    ):
        if len(files) < 1:
            raise ValueError("At least one file is required")
        mode = "a" if append else "w"
        self.files = files
        self.number_of_files = len(files)
        self.writers: List[Union[FastaWriter, FastqWriter]] = [
            _open_single(  # type: ignore
                file,
                opener=opener,
                fileformat="fasta",
                mode=mode,
                qualities=False,
            )
            for file in self.files
        ]

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}"
            f"({', '.join(repr(writer) for writer in self.writers)})"
        )

    def close(self):
        for writer in self.writers:
            writer.close()

    def write(self, *records: SequenceRecord):
        if len(records) != self.number_of_files:
            raise ValueError(f"records must have length {self.number_of_files}")
        for record, writer in zip(records, self.writers):
            writer.write(record)

    def write_iterable(self, records_iterable: Iterable[Tuple[SequenceRecord, ...]]):
        for records in records_iterable:
            self.write(*records)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class MultipleFastqWriter(MultipleFileWriter):
    def __init__(
        self,
        *files: Union[str, PathLike, BinaryIO],
        opener=xopen,
        append: bool = False,
    ):
        if len(files) < 1:
            raise ValueError("At least one file is required")
        mode = "a" if append else "w"
        self.files = files
        self.number_of_files = len(files)
        self.writers: List[IO] = [
            opener(file, mode + "b") if not hasattr(file, "write") else file
            for file in self.files
        ]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}" f"({', '.join(str(f) for f in self.files)})"

    def close(self):
        for writer in self.writers:
            writer.close()

    def write(self, *records: SequenceRecord):
        if len(records) != self.number_of_files:
            raise ValueError(f"records must have length {self.number_of_files}")
        for record, writer in zip(records, self.writers):
            writer.write(record.fastq_bytes())

    def write_iterable(self, records_iterable: Iterable[Tuple[SequenceRecord, ...]]):
        # Use faster methods for more common cases before falling back to
        # generic multiple files mode (which is much slower due to calling the
        # zip function).
        if self.number_of_files == 1:
            output = self.writers[0]
            for (record,) in records_iterable:
                output.write(record.fastq_bytes())
        elif self.number_of_files == 2:
            output1 = self.writers[0]
            output2 = self.writers[1]
            for record1, record2 in records_iterable:
                output1.write(record1.fastq_bytes())
                output2.write(record2.fastq_bytes())
        elif self.number_of_files == 3:
            output1 = self.writers[0]
            output2 = self.writers[1]
            output3 = self.writers[2]
            for record1, record2, record3 in records_iterable:
                output1.write(record1.fastq_bytes())
                output2.write(record2.fastq_bytes())
                output3.write(record3.fastq_bytes())
        else:  # More than 3 files is quite uncommon.
            writers = self.writers
            for records in records_iterable:
                for record, output in zip(records, writers):
                    output.write(record.fastq_bytes())

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
