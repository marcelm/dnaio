from abc import ABC, abstractmethod
from typing import Iterator, Tuple

from dnaio import SequenceRecord


class SingleEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[SequenceRecord]:
        """
        Iterate over an input containing sequence records

        Yields:
            `SequenceRecord` objects

        Raises:
            `FileFormatError` if there was a parse error
        """


class PairedEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
        """
        Iterate over an input containing paired-end records

        Yields:
            Pairs of `SequenceRecord` objects

        Raises:
            `FileFormatError` if there was a parse error or if reads are improperly paired,
            that is, if there are more reads in one file than the other or if the record IDs
            do not match (according to `SequenceRecord.is_mate`).
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
        this method.
        """
