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


def test_write(tmpdir):
    s = dnaio.Sequence('name', 'ACGT', 'HHHH')
    out_fastq = tmpdir.join('out.fastq')
    with dnaio.open(str(out_fastq), mode='w') as f:
        f.write(s)
    assert out_fastq.read() == '@name\nACGT\n+\nHHHH\n'


def test_write_gz(tmpdir):
    s = dnaio.Sequence('name', 'ACGT', 'HHHH')
    out_fastq = tmpdir.join('out.fastq.gz')
    with dnaio.open(str(out_fastq), mode='w') as f:
        f.write(s)

    import gzip
    with gzip.open(str(out_fastq)) as f:
        assert f.read() == b'@name\nACGT\n+\nHHHH\n'


def test_write_gz_with_xopen(tmpdir):
    s = dnaio.Sequence('name', 'ACGT', 'HHHH')
    out_fastq = tmpdir.join('out.fastq.gz')
    with xopen(str(out_fastq), 'wb') as gzf:
        with dnaio.open(gzf, mode='w') as f:
            f.write(s)

    import gzip
    with gzip.open(str(out_fastq)) as f:
        assert f.read() == b'@name\nACGT\n+\nHHHH\n'


def test_read_pathlib_path():
    path = Path('tests/data/simple.fasta')
    with dnaio.open(path) as f:
        records = list(f)
        assert len(records) == 2
