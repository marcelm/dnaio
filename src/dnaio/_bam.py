import struct
from io import BufferedIOBase


def read_bam_header(fileobj: BufferedIOBase) -> bytes:
    magic = fileobj.read(4)
    if not isinstance(magic, bytes):
        raise TypeError(
            f"fileobj {fileobj} is not a binary IO type, " f"got {type(fileobj)}"
        )
    if len(magic) < 4:
        raise EOFError("Truncated BAM file")
    if magic[:4] != b"BAM\1":
        raise ValueError(
            f"fileobj: {fileobj}, is not a BAM file. No BAM magic, instead "
            f"found {magic!r}"
        )
    return read_bam_header_after_magic(fileobj)


def read_bam_header_after_magic(fileobj: BufferedIOBase) -> bytes:
    header_size = fileobj.read(4)
    if len(header_size) < 4:
        raise EOFError("Truncated BAM file")
    (l_text,) = struct.unpack("<I", header_size)
    header = fileobj.read(l_text)
    if len(header) < l_text:
        raise EOFError("Truncated BAM file")
    n_ref_obj = fileobj.read(4)
    if len(n_ref_obj) < 4:
        raise EOFError("Truncated BAM file")
    (n_ref,) = struct.unpack("<I", n_ref_obj)
    for i in range(n_ref):
        l_name_obj = fileobj.read(4)
        if len(l_name_obj) < 4:
            raise EOFError("Truncated BAM file")
        (l_name,) = struct.unpack("<I", l_name_obj)
        reference_chunk_size = l_name + 4  # Include name and uint32_t of size
        reference_chunk = fileobj.read(reference_chunk_size)
        if len(reference_chunk) < reference_chunk_size:
            raise EOFError("Truncated BAM file")
    # Fileobj is now skipped ahead and at the start of the BAM records
    return header
