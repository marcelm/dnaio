from abc import ABC, abstractmethod
from typing import Iterable, Iterator, Tuple

from dnaio import SequenceRecord


class SingleEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[SequenceRecord]:
        """Yield the records in the input as `SequenceRecord` objects."""


class PairedEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
        """
        Yield the records in the paired-end input as pairs of `SequenceRecord` objects.

        Raises a `FileFormatError` if reads are improperly paired, that is,
        if there are more reads in one file than the other or if the record IDs
        do not match (according to `SequenceRecord.is_mate`).
        """


class MultipleFileReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[SequenceRecord, ...]]:
        """
        Yields the records as N-tuples depending on the number of input files.

        If there are two or more files, it will raise a `FileFormatError`
        if reads are improperly paired. That is when there are more reads in
        one file than the others or if the record IDs do not match (according
        to `SequenceRecord.is_mate` or `records_are_mates`).
        """


class SingleEndWriter(ABC):
    @abstractmethod
    def write(self, record: SequenceRecord) -> None:
        """Write a `SequenceRecord` to the output."""


class PairedEndWriter(ABC):
    @abstractmethod
    def write(self, record1: SequenceRecord, record2: SequenceRecord) -> None:
        """
        Write a pair of `SequenceRecord` objects to the paired-end output.

        This method does not verify that both records have matching IDs
        because this was already done at parsing time. If it is possible
        that the record IDs no longer match, check that
        ``record1.is_mate(record2)`` returns True before calling
        this function.
        """


class MultipleFileWriter(ABC):
    number_of_files: int

    @abstractmethod
    def write(self, records: Tuple[SequenceRecord, ...]) -> None:
        """
        Write a N-tuple of SequenceRecords to the output. N should be equal
        to the number of files the MultipleFileWriter was initiated with.

        This method does not check wheter the records are properly paired.
        """

    def write_list(self, list_of_records: Iterable[Tuple[SequenceRecord, ...]]):
        """
        Iterates over the list (or other iterable container) and writes all
        N-tuples of SequenceRecord to disk. N should be equal
        to the number of files the MultipleFileWriter was initiated with.

        This method does not check wheter the records are properly paired.
        """
