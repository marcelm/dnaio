import sys
from typing import (
    Generic,
    Optional,
    Tuple,
    BinaryIO,
    Iterator,
    Type,
    TypeVar,
    ByteString,
)

class SequenceRecord:
    name: str
    sequence: str
    qualities: Optional[str]
    def __init__(
        self, name: str, sequence: str, qualities: Optional[str] = ...
    ) -> None: ...
    def __getitem__(self, s: slice) -> SequenceRecord: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...
    def __richcmp__(self, other: SequenceRecord, op: int) -> bool: ...
    def qualities_as_bytes(self) -> bytes: ...
    def fastq_bytes(self, two_headers: bool = ...) -> bytes: ...
    def is_mate(self, other: SequenceRecord) -> bool: ...
    def reverse_complement(self) -> SequenceRecord: ...
    @property
    def id(self) -> str: ...
    @property
    def comment(self) -> Optional[str]: ...

# Bytestring = Union[bytes, bytearray, memoryview]. Technically just 'bytes' is
# acceptable as an alias, but even more technically this function supports all
# types that implement the buffer protocol, for which there is no type yet.
# See: https://github.com/python/typing/issues/593
def paired_fastq_heads(
    buf1: ByteString, buf2: ByteString, end1: int, end2: int
) -> Tuple[int, int]: ...
def bam_head(buf: bytes, end: Optional[int] = None) -> int: ...
def records_are_mates(
    __first_record: SequenceRecord,
    __second_record: SequenceRecord,
    *__other_records: SequenceRecord,
) -> bool: ...

T = TypeVar("T")

class FastqIter(Generic[T]):
    def __init__(
        self, file: BinaryIO, sequence_class: Type[T], buffer_size: int = ...
    ): ...
    def __iter__(self) -> Iterator[T]: ...
    def __next__(self) -> T: ...
    @property
    def number_of_records(self) -> int: ...

class BamIter:
    def __init__(self, file: BinaryIO, read_in_size: int, with_header: bool = True): ...
    def __iter__(self) -> Iterator[SequenceRecord]: ...
    def __next__(self) -> SequenceRecord: ...
    @property
    def header(self) -> bytes: ...
    @property
    def number_of_records(self) -> int: ...

# Deprecated
def record_names_match(header1: str, header2: str) -> bool: ...

# Private
def bytes_ascii_check(b: bytes, length: int = -1) -> bool: ...
