import itertools
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
