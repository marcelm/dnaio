from os import PathLike
from typing import BinaryIO, Iterable, Iterator, List, Optional, Tuple, Union

from xopen import xopen

from ._core import SequenceRecord, records_are_mates
from .exceptions import FileFormatError
from .readers import FastaReader, FastqReader
from .singleend import _open_single
from .writers import FastaWriter, FastqWriter


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
            _open_single(
                file, opener=opener, fileformat=fileformat, mode="r"
            )  # type: ignore
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
            return zip(self.readers[0])
        elif len(self.files) == 2:
            for record1, record2 in zip(*self.readers):
                if not record1.is_mate(record2):
                    raise FileFormatError(
                        f"Records are out of sync, names "
                        f"{record1}, f{record2} do not match.",
                        line=None,
                    )
                yield record1, record2
        else:
            for records in zip(*self.readers):
                if not records_are_mates(*records):
                    raise FileFormatError(
                        f"Records are out of sync, names "
                        f"{', '.join(r.name for r in records)} do not match.",
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


class MultipleFileWriter:
    def __init__(
        self,
        *files: Union[str, PathLike, BinaryIO],
        fileformat: Optional[str] = "fastq",
        qualities: Optional[bool] = None,
        opener=xopen,
        append: bool = False,
    ):
        if len(files) < 1:
            raise ValueError("At least one file is required")
        mode = "a" if append else "w"
        self.files = files
        self.number_of_files = len(files)
        self.writers: List[Union[FastaWriter, FastqWriter]] = [
            _open_single(
                file,
                opener=opener,
                fileformat=fileformat,
                mode=mode,
                qualities=qualities,
            )  # type: ignore
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

    def write(self, records: Tuple[SequenceRecord, ...]):
        if len(records) != self.number_of_files:
            raise ValueError(f"records must have length {self.number_of_files}")
        for record, writer in zip(records, self.writers):
            writer.write(record)

    def write_iterable(self,
                       records_iterable: Iterable[Tuple[SequenceRecord, ...]]):
        for records in records_iterable:
            self.write(records)
