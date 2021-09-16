from abc import ABC, abstractmethod
from typing import Iterator, Tuple

from dnaio import Sequence


class SingleEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Sequence]:
        pass


class PairedEndReader(ABC):
    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[Sequence, Sequence]]:
        pass


class SingleEndWriter(ABC):
    @abstractmethod
    def write(self, record: Sequence) -> None:
        pass


class PairedEndWriter(ABC):
    @abstractmethod
    def write(self, record1: Sequence, record2: Sequence) -> None:
        pass
