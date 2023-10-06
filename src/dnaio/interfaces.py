from abc import abstractmethod
from contextlib import AbstractContextManager
from typing import Iterable, Iterator, Tuple

from dnaio import SequenceRecord


class SingleEndReader(AbstractContextManager):
    delivers_qualities: bool
    number_of_records: int

    @abstractmethod
    def __iter__(self) -> Iterator[SequenceRecord]:
        """
        Iterate over an input containing sequence records

        Yields:
            `SequenceRecord` objects

        Raises:
            `FileFormatError`
                if there was a parse error
        """

    @abstractmethod
    def close(self) -> None:
        pass


class PairedEndReader(AbstractContextManager):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
        """
        Iterate over an input containing paired-end records

        Yields:
            Pairs of `SequenceRecord` objects

        Raises:
            `FileFormatError`
                if there was a parse error or if reads are improperly paired,
                that is, if there are more reads in one file than the other or
                if the record IDs do not match (according to
                `SequenceRecord.is_mate`).
        """

    @abstractmethod
    def close(self) -> None:
        pass


class SingleEndWriter(AbstractContextManager):
    @abstractmethod
    def write(self, record: SequenceRecord) -> None:
        """Write a `SequenceRecord` to the output."""

    @abstractmethod
    def close(self) -> None:
        pass


class PairedEndWriter(AbstractContextManager):
    @abstractmethod
    def write(self, record1: SequenceRecord, record2: SequenceRecord) -> None:
        """
        Write a pair of `SequenceRecord` objects to the paired-end output.

        This method does not verify that both records have matching IDs
        because this was already done at parsing time. If it is possible
        that the record IDs no longer match, check that
        ``record1.is_mate(record2)`` returns True before calling
        this method.
        """

    @abstractmethod
    def close(self) -> None:
        pass


class MultipleFileWriter(AbstractContextManager):
    _number_of_files: int

    @abstractmethod
    def write(self, *records: SequenceRecord) -> None:
        """
        Write N SequenceRecords to the output. N must be equal
        to the number of files the MultipleFileWriter was initialized with.

        This method does not check whether the records are properly paired.
        """

    @abstractmethod
    def write_iterable(self, list_of_records: Iterable[Tuple[SequenceRecord, ...]]):
        """
        Iterate over the list (or other iterable container) and write all
        N-tuples of SequenceRecord to disk. N must be equal
        to the number of files the MultipleFileWriter was initialized with.

        This method does not check whether the records are properly paired.
        This method may provide a speed boost over calling write for each
        tuple of SequenceRecords individually.
        """

    @abstractmethod
    def close(self) -> None:
        pass
