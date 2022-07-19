import io
import itertools
import os
from pathlib import Path

from dnaio import SequenceRecord, open_multiple

import pytest


@pytest.mark.parametrize(
    ["fileformat", "number_of_files"],
    itertools.product(("fasta", "fastq"), (1, 2, 3, 4)),
)
def test_read_files(fileformat, number_of_files):
    file = Path(__file__).parent / "data" / ("simple." + fileformat)
    files = [file for _ in range(number_of_files)]
    with open_multiple(*files) as multiple_reader:
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
        open_multiple(**kwargs)


def test_open_multiple_unsupported_mode():
    with pytest.raises(ValueError) as error:
        open_multiple(os.devnull, mode="X")
    error.match("mode")


@pytest.mark.parametrize(
    ["number_of_files", "content"],
    itertools.product(
        (1, 2, 3, 4), (">my_fasta\nAGCTAGA\n", "@my_fastq\nAGC\n+\nHHH\n")
    ),
)
def test_multiple_binary_read(number_of_files, content):
    files = [io.BytesIO(content.encode("ascii")) for _ in range(number_of_files)]
    with open_multiple(*files) as reader:
        for records_tup in reader:
            pass


@pytest.mark.parametrize(
    ["number_of_files", "fileformat"],
    itertools.product((1, 2, 3, 4), ("fastq", "fasta")),
)
def test_multiple_binary_write(number_of_files, fileformat):
    files = [io.BytesIO() for _ in range(number_of_files)]
    records = [SequenceRecord("A", "A", "A") for _ in range(number_of_files)]
    with open_multiple(*files, mode="w", fileformat=fileformat) as writer:
        writer.write(*records)


@pytest.mark.parametrize(
    ["number_of_files", "fileformat"],
    itertools.product((1, 2, 3, 4), ("fastq", "fasta")),
)
def test_multiple_write_too_much(number_of_files, fileformat):
    files = [io.BytesIO() for _ in range(number_of_files)]
    records = [SequenceRecord("A", "A", "A") for _ in range(number_of_files + 1)]
    with open_multiple(*files, mode="w", fileformat=fileformat) as writer:
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
    with open_multiple(*files, mode="w", fileformat=fileformat) as writer:
        writer.write_iterable(records_list)
