from typing import Optional, Tuple, BinaryIO, Iterator, Type, TypeVar, ByteString


class Sequence:
    name: str
    sequence: str
    qualities: Optional[str]
    def __init__(self, name: str, sequence: str, qualities: Optional[str] = ...) -> None: ...
    def __getitem__(self, s: slice) -> Sequence: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...
    def __richcmp__(self, other: Sequence, op: int) -> bool: ...
    def qualities_as_bytes(self) -> bytes: ...
    def fastq_bytes(self) -> bytes: ...
    def fastq_bytes_two_headers(self) -> bytes: ...
    def is_mate(self, other: Sequence) -> bool: ...

class BytesSequence:
    name: bytes
    sequence: bytes
    qualities: bytes

    def __init__(self, name: bytes, sequence: bytes, qualities: bytes = ...) -> None: ...
    def __getitem__(self, s: slice) -> BytesSequence: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...
    def __richcmp__(self, other: BytesSequence, op: int) -> bool: ...
    def fastq_bytes(self) -> bytes: ...
    def fastq_bytes_two_headers(self) -> bytes: ...
    def is_mate(self, other: BytesSequence) -> bool: ...

# Bytestring = Union[bytes, bytearray, memoryview]. Technically just 'bytes' is
# acceptable as an alias, but even more technically this function supports all
# types that implement the buffer protocol, for which there is no type yet.
# See: https://github.com/python/typing/issues/593
def paired_fastq_heads(buf1: ByteString, buf2: ByteString, end1: int, end2: int) -> Tuple[int, int]: ...
def record_names_match(header1: str, header2: str) -> bool: ...
def record_names_match_bytes(header1: bytes, header2: bytes) -> bool: ...

T = TypeVar("T")
class FastqIter:
    def __init__(self, file: BinaryIO, sequence_class: Type[T], buffer_size: int = ...): ...
    def __iter__(self) -> Iterator[T]: ...
    def __next__(self) -> T: ...

    @property
    def n_records(self) -> int: ...
