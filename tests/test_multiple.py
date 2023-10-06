import io
import itertools
import os
from pathlib import Path

import dnaio
from dnaio import SequenceRecord, _open_multiple

import pytest


@pytest.mark.parametrize(
    ["fileformat", "number_of_files"],
    itertools.product(("fasta", "fastq"), (1, 2, 3, 4)),
)
def test_read_files(fileformat, number_of_files):
    file = Path(__file__).parent / "data" / ("simple." + fileformat)
    files = [file] * number_of_files
    with _open_multiple(*files) as multiple_reader:
        for records in multiple_reader:
            pass
        assert len(records) == number_of_files
        assert isinstance(records, tuple)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(mode="w", fileformat="fasta"),
        dict(mode="r"),
        dict(mode="w", fileformat="fastq"),
    ],
)
def test_open_no_file_error(kwargs):
    with pytest.raises(ValueError):
        _open_multiple(**kwargs)


def test_open_multiple_unsupported_mode():
    with pytest.raises(ValueError) as error:
        _open_multiple(os.devnull, mode="X")
    error.match("one of 'r', 'w', 'a'")


@pytest.mark.parametrize(
    ["number_of_files", "content"],
    itertools.product(
        (1, 2, 3, 4), (">my_fasta\nAGCTAGA\n", "@my_fastq\nAGC\n+\nHHH\n")
    ),
)
def test_multiple_binary_read(number_of_files, content):
    files = [io.BytesIO(content.encode("ascii")) for _ in range(number_of_files)]
    with _open_multiple(*files) as reader:
        for records_tup in reader:
            pass


@pytest.mark.parametrize(
    ["number_of_files", "fileformat"],
    itertools.product((1, 2, 3, 4), ("fastq", "fasta")),
)
def test_multiple_binary_write(number_of_files, fileformat):
    files = [io.BytesIO() for _ in range(number_of_files)]
    records = [SequenceRecord("A", "A", "A") for _ in range(number_of_files)]
    with _open_multiple(*files, mode="w", fileformat=fileformat) as writer:
        writer.write(*records)


@pytest.mark.parametrize(
    ["number_of_files", "fileformat"],
    itertools.product((1, 2, 3, 4), ("fastq", "fasta")),
)
def test_multiple_write_too_many(number_of_files, fileformat):
    files = [io.BytesIO() for _ in range(number_of_files)]
    records = [SequenceRecord("A", "A", "A") for _ in range(number_of_files + 1)]
    with _open_multiple(*files, mode="w", fileformat=fileformat) as writer:
        with pytest.raises(ValueError) as error:
            writer.write(*records)
    error.match(str(number_of_files))


@pytest.mark.parametrize(
    ["number_of_files", "fileformat"],
    itertools.product((1, 2, 3, 4), ("fastq", "fasta")),
)
def test_multiple_write_iterable(number_of_files, fileformat):
    files = [io.BytesIO() for _ in range(number_of_files)]
    records = [SequenceRecord("A", "A", "A") for _ in range(number_of_files)]
    records_list = [records, records, records]
    with _open_multiple(*files, mode="w", fileformat=fileformat) as writer:
        writer.write_iterable(records_list)


@pytest.mark.parametrize("number_of_files", (2, 3, 4))
def test_multiple_read_unmatched_names(number_of_files):
    record1_content = b"@my_fastq\nAGC\n+\nHHH\n"
    record2_content = b"@my_fasterq\nAGC\n+\nHHH\n"
    files = (
        io.BytesIO(record1_content),
        *(io.BytesIO(record2_content) for _ in range(number_of_files - 1)),
    )
    with _open_multiple(*files) as reader:
        with pytest.raises(dnaio.FileFormatError) as error:
            for records in reader:
                pass  # pragma: no coverage
    error.match("do not match")


@pytest.mark.parametrize("number_of_files", (2, 3, 4))
def test_multiple_read_out_of_sync(number_of_files):
    record1_content = b"@my_fastq\nAGC\n+\nHHH\n"
    record2_content = b"@my_fastq\nAGC\n+\nHHH\n@my_secondfastq\nAGC\n+\nHHH\n"
    files = (
        io.BytesIO(record1_content),
        *(io.BytesIO(record2_content) for _ in range(number_of_files - 1)),
    )
    with _open_multiple(*files) as reader:
        with pytest.raises(dnaio.FileFormatError) as error:
            for records in reader:
                pass
    error.match("unequal amount")
