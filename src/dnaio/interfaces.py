from abc import ABC, abstractmethod
from typing import Iterator, Tuple

from dnaio import SequenceRecord


class SingleEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[SequenceRecord]:
        pass


class PairedEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[SequenceRecord, SequenceRecord]]:
        pass


class SingleEndWriter(ABC):
    @abstractmethod
    def write(self, record: SequenceRecord) -> None:
        pass


class PairedEndWriter(ABC):
    @abstractmethod
    def write(self, record1: SequenceRecord, record2: SequenceRecord) -> None:
        pass
