import itertools
import os
from pathlib import Path
from dnaio.multipleend import open_multiple

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

