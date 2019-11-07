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


def formatted_sequence(record, fileformat):
    if fileformat == "fastq":
        return "@{}\n{}\n+\n{}\n".format(record.name, record.sequence, record.qualities)
    else:
        return ">{}\n{}\n".format(record.name, record.sequence)


def formatted_sequences(records, fileformat):
    return "".join(formatted_sequence(record, fileformat) for record in records)


def test_formatted_sequence():
    s = dnaio.Sequence("s1", "ACGT", "HHHH")
    assert ">s1\nACGT\n" == formatted_sequence(s, "fasta")
    assert "@s1\nACGT\n+\nHHHH\n" == formatted_sequence(s, "fastq")


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


def test_read_opener(fileformat, extension):
    def my_opener(path, mode):
        import io
        if fileformat == "fasta":
            data = b">read\nACG\n"
        else:
            data = b"@read\nACG\n+\nHHH\n"
        return io.BytesIO(data)

    with dnaio.open("totally-ignored-filename." + fileformat + extension, opener=my_opener) as f:
        records = list(f)
    assert len(records) == 1
    assert records[0].name == "read"
    assert records[0].sequence == "ACG"


@pytest.mark.parametrize("interleaved", [False, True])
def test_paired_opener(fileformat, extension, interleaved):
    def my_opener(_path, _mode):
        import io
        if fileformat == "fasta":
            data = b">read\nACG\n"
        else:
            data = b"@read\nACG\n+\nHHH\n"
        return io.BytesIO(data + data)

    path1 = "ignored-filename." + fileformat + extension
    path2 = "also-ignored-filename." + fileformat + extension
    if interleaved:
        with dnaio.open(path1, file2=path2, opener=my_opener) as f:
            records = list(f)
        expected = 2
    else:
        with dnaio.open(path1, interleaved=True, opener=my_opener) as f:
            records = list(f)
        expected = 1
    assert len(records) == expected
    assert records[0][0].name == "read"
    assert records[0][0].sequence == "ACG"
    assert records[0][1].name == "read"
    assert records[0][1].sequence == "ACG"


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
    path = Path(str(tmpdir / ("out." + fileformat + extension)))
    with dnaio.open(path, mode="w") as f:
        f.write(s1)
    if fileformat == "fasta":
        expected = b">s1\nACGT\n"
    else:
        expected = b"@s1\nACGT\n+\nHHHH\n"
    with xopen(path, "rb") as f:
        assert f.read() == expected


def test_write_paired_same_path(tmpdir):
    path1 = str(tmpdir / "same.fastq")
    path2 = str(tmpdir / "same.fastq")
    with pytest.raises(ValueError):
        with dnaio.open(file1=path1, file2=path2, mode="w"):
            pass


def test_write_paired(tmpdir, fileformat, extension):
    r1 = [
        dnaio.Sequence("s1", "ACGT", "HHHH"),
        dnaio.Sequence("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.Sequence("t1", "TCGT", "5HHH"),
        dnaio.Sequence("t2", "TGCA", "5383"),
    ]
    path1 = str(tmpdir / ("out.1." + fileformat + extension))
    path2 = str(tmpdir / ("out.2." + fileformat + extension))

    with dnaio.open(path1, file2=path2, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
        f.write(r1[1], r2[1])
    with xopen(path1) as f:
        assert formatted_sequences(r1, fileformat) == f.read()
    with xopen(path2) as f:
        assert formatted_sequences(r2, fileformat) == f.read()


def test_write_interleaved(tmpdir, fileformat, extension):
    r1 = [
        dnaio.Sequence("s1", "ACGT", "HHHH"),
        dnaio.Sequence("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.Sequence("t1", "TCGT", "5HHH"),
        dnaio.Sequence("t2", "TGCA", "5383"),
    ]
    path = str(tmpdir / ("out.interleaved." + fileformat + extension))

    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
        f.write(r1[1], r2[1])
    expected = [r1[0], r2[0], r1[1], r2[1]]
    with xopen(path) as f:
        assert formatted_sequences(expected, fileformat) == f.read()


def test_append(tmpdir, fileformat, extension):
    s1 = dnaio.Sequence("s1", "ACGT", "HHHH")
    s2 = dnaio.Sequence("s2", "CGCA", "8383")
    path = str(tmpdir / ("out." + fileformat + extension))
    with dnaio.open(path, mode="w") as f:
        f.write(s1)
    with dnaio.open(path, mode="a") as f:
        f.write(s2)
    with xopen(path) as f:
        assert formatted_sequences([s1, s2], fileformat) == f.read()


def test_append_paired(tmpdir, fileformat, extension):
    r1 = [
        dnaio.Sequence("s1", "ACGT", "HHHH"),
        dnaio.Sequence("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.Sequence("t1", "TCGT", "5HHH"),
        dnaio.Sequence("t2", "TGCA", "5383"),
    ]
    path1 = str(tmpdir / ("out.1." + fileformat + extension))
    path2 = str(tmpdir / ("out.2." + fileformat + extension))

    with dnaio.open(path1, file2=path2, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
    with dnaio.open(path1, file2=path2, fileformat=fileformat, mode="a") as f:
        f.write(r1[1], r2[1])
    with xopen(path1) as f:
        assert formatted_sequences(r1, fileformat) == f.read()
    with xopen(path2) as f:
        assert formatted_sequences(r2, fileformat) == f.read()


def test_append_interleaved(tmpdir, fileformat, extension):
    r1 = [
        dnaio.Sequence("s1", "ACGT", "HHHH"),
        dnaio.Sequence("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.Sequence("t1", "TCGT", "5HHH"),
        dnaio.Sequence("t2", "TGCA", "5383"),
    ]
    path = str(tmpdir / ("out.interleaved." + fileformat + extension))

    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="a") as f:
        f.write(r1[1], r2[1])
    expected = [r1[0], r2[0], r1[1], r2[1]]
    with xopen(path) as f:
        assert formatted_sequences(expected, fileformat) == f.read()


def make_random_fasta(path, n_records):
    from random import choice
    with xopen(path, "w") as f:
        for i in range(n_records):
            name = "sequence_{}".format(i)
            sequence = "".join(choice("ACGT") for _ in range(300))
            print(">", name, "\n", sequence, sep="", file=f)


def test_islice_gzip_does_not_fail(tmpdir):
    path = str(tmpdir / "file.fasta.gz")
    make_random_fasta(path, 100)
    f = dnaio.open(path)
    next(iter(f))
    f.close()
