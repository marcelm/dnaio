from pathlib import Path

import dnaio
from xopen import xopen

import pytest


@pytest.fixture(params=["", ".gz", ".bz2", ".xz"])
def extension(request):
    return request.param


@pytest.fixture(params=["fasta", "fastq"])
def fileformat(request):
    return request.param


SIMPLE_RECORDS = {
    "fasta": [
        dnaio.Sequence("first_sequence", "SEQUENCE1"),
        dnaio.Sequence("second_sequence", "SEQUENCE2"),
    ],
    "fastq": [
        dnaio.Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
        dnaio.Sequence("second_sequence", "SEQUENCE2", "83<??:(61"),
    ],
}


def test_version():
    _ = dnaio.__version__


def test_read(fileformat, extension):
    with dnaio.open("tests/data/simple." + fileformat + extension) as f:
        records = list(f)
    assert records == SIMPLE_RECORDS[fileformat]


def test_read_pathlib_path(fileformat, extension):
    path = Path("tests/data/simple." + fileformat + extension)
    with dnaio.open(path) as f:
        records = list(f)
    assert records == SIMPLE_RECORDS[fileformat]


def test_detect_fastq_from_content():
    """FASTQ file that is not named .fastq"""
    with dnaio.open('tests/data/missingextension') as f:
        record = next(iter(f))
        assert record.name == 'prefix:1_13_573/1'


def test_detect_compressed_fastq_from_content():
    """Compressed FASTQ file that is not named .fastq.gz"""
    with dnaio.open('tests/data/missingextension.gz') as f:
        record = next(iter(f))
    assert record.name == 'prefix:1_13_573/1'


def test_write(tmpdir, extension):
    s = dnaio.Sequence('name', 'ACGT', 'HHHH')
    out_fastq = tmpdir.join("out.fastq" + extension)
    with dnaio.open(str(out_fastq), mode='w') as f:
        f.write(s)
    with xopen(out_fastq) as f:
        assert f.read() == '@name\nACGT\n+\nHHHH\n'


def test_write_with_xopen(tmpdir, fileformat, extension):
    s = dnaio.Sequence('name', 'ACGT', 'HHHH')
    out_fastq = str(tmpdir.join("out." + fileformat + extension))
    with xopen(out_fastq, 'wb') as outer_f:
        with dnaio.open(outer_f, mode='w', fileformat=fileformat) as f:
            f.write(s)

    with xopen(out_fastq) as f:
        if fileformat == "fasta":
            assert f.read() == ">name\nACGT\n"
        else:
            assert f.read() == "@name\nACGT\n+\nHHHH\n"


def test_write_pathlib(tmpdir, fileformat, extension):
    s1 = dnaio.Sequence("s1", "ACGT", "HHHH")
    path = Path(tmpdir / ("out." + fileformat + extension))
    with dnaio.open(path, mode="w") as f:
        f.write(s1)
    if fileformat == "fasta":
        expected = b">s1\nACGT\n"
    else:
        expected = b"@s1\nACGT\n+\nHHHH\n"
    with xopen(path, "rb") as f:
        assert f.read() == expected
