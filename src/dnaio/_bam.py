from typing import BinaryIO


def read_bam_header(fileobj: BinaryIO) -> bytes:
    magic_and_header_size = fileobj.read(8)
    if not isinstance(magic_and_header_size, bytes):
        raise TypeError(
            f"fileobj {fileobj} is not a binary IO type, " f"got {type(fileobj)}"
        )
    if len(magic_and_header_size) < 8:
        raise EOFError("Truncated BAM file")
    if magic_and_header_size[:4] != b"BAM\1":
        raise ValueError(
            f"fileobj: {fileobj}, is not a BAM file. No BAM magic, instead "
            f"found {magic_and_header_size[:4]}"
        )
    l_text = int.from_bytes(magic_and_header_size[4:], "little", signed=False)
    header = fileobj.read(l_text)
    if len(header) < l_text:
        raise EOFError("Truncated BAM file")
    n_ref_obj = fileobj.read(4)
    if len(n_ref_obj) < 4:
        raise EOFError("Truncated BAM file")
    n_ref = int.from_bytes(n_ref_obj, "little", signed=False)
    for i in range(n_ref):
        l_name_obj = fileobj.read(4)
        if len(l_name_obj) < 4:
            raise EOFError("Truncated BAM file")
        l_name = int.from_bytes(l_name_obj, "little", signed=False)
        reference_chunk_size = l_name + 4  # Include name and uint32_t of size
        reference_chunk = fileobj.read(reference_chunk_size)
        if len(reference_chunk) < reference_chunk_size:
            raise EOFError("Truncated BAM file")
    # Fileobj is now skipped ahead and at the start of the BAM records
    return header
