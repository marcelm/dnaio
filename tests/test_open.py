import os
from pathlib import Path

import pytest
from xopen import xopen

import dnaio
from dnaio import FileFormatError, UnknownFileFormat


@pytest.fixture(params=["", ".gz", ".bz2", ".xz"])
def extension(request):
    return request.param


@pytest.fixture(params=["fasta", "fastq"])
def fileformat(request):
    return request.param


SIMPLE_RECORDS = {
    "fasta": [
        dnaio.SequenceRecord("first_sequence", "SEQUENCE1"),
        dnaio.SequenceRecord("second_sequence", "SEQUENCE2"),
    ],
    "fastq": [
        dnaio.SequenceRecord("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
        dnaio.SequenceRecord("second_sequence", "SEQUENCE2", "83<??:(61"),
    ],
}


def formatted_sequence(record, fileformat) -> str:
    if fileformat == "fastq":
        return "@{}\n{}\n+\n{}\n".format(record.name, record.sequence, record.qualities)
    else:
        return ">{}\n{}\n".format(record.name, record.sequence)


def formatted_sequences(records, fileformat) -> str:
    return "".join(formatted_sequence(record, fileformat) for record in records)


def test_formatted_sequence() -> None:
    s = dnaio.SequenceRecord("s1", "ACGT", "HHHH")
    assert ">s1\nACGT\n" == formatted_sequence(s, "fasta")
    assert "@s1\nACGT\n+\nHHHH\n" == formatted_sequence(s, "fastq")


def test_version() -> None:
    _ = dnaio.__version__


def test_open_nonexistent(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        with dnaio.open(tmp_path / "nonexistent"):
            pass  # pragma: no cover


def test_open_empty_file_with_unrecognized_extension(tmp_path) -> None:
    path = tmp_path / "unrecognized-extension.tmp"
    path.touch()
    with dnaio.open(path) as f:
        records = list(f)
    assert records == []


def test_fileformat_error(tmp_path) -> None:
    with open(tmp_path / "file.fastq", mode="w") as f:
        print("this is not a FASTQ file", file=f)
    with pytest.raises(FileFormatError) as e:
        with dnaio.open(tmp_path / "file.fastq") as f:
            _ = list(f)  # pragma: no cover
    assert "at line 2" in str(e.value)  # Premature end of file


def test_write_unknown_file_format(tmp_path) -> None:
    with pytest.raises(UnknownFileFormat):
        with dnaio.open(tmp_path / "out.txt", mode="w") as f:
            f.write(dnaio.SequenceRecord("name", "ACG", "###"))  # pragma: no cover


def test_read_unknown_file_format(tmp_path) -> None:
    with open(tmp_path / "file.txt", mode="w") as f:
        print("text file", file=f)
    with pytest.raises(UnknownFileFormat):
        with dnaio.open(tmp_path / "file.txt") as f:
            _ = list(f)  # pragma: no cover


def test_invalid_format(tmp_path) -> None:
    with pytest.raises(UnknownFileFormat):
        with dnaio.open(tmp_path / "out.txt", mode="w", fileformat="foo"):
            pass  # pragma: no cover


def test_write_qualities_to_file_without_fastq_extension(tmp_path) -> None:
    with dnaio.open(tmp_path / "out.txt", mode="w", qualities=True) as f:
        f.write(dnaio.SequenceRecord("name", "ACG", "###"))

    with dnaio.open(tmp_path / "out.txt", mode="w", qualities=False) as f:
        f.write(dnaio.SequenceRecord("name", "ACG", None))


def test_read(fileformat, extension) -> None:
    with dnaio.open("tests/data/simple." + fileformat + extension) as f:
        records = list(f)
    assert records == SIMPLE_RECORDS[fileformat]


def test_read_pathlib_path(fileformat, extension) -> None:
    path = Path("tests/data/simple." + fileformat + extension)
    with dnaio.open(path) as f:
        records = list(f)
    assert records == SIMPLE_RECORDS[fileformat]


def test_read_opener(fileformat, extension) -> None:
    def my_opener(path, mode):
        import io

        if fileformat == "fasta":
            data = b">read\nACG\n"
        else:
            data = b"@read\nACG\n+\nHHH\n"
        return io.BytesIO(data)

    with dnaio.open(
        "totally-ignored-filename." + fileformat + extension, opener=my_opener
    ) as f:
        records = list(f)
    assert len(records) == 1
    assert records[0].name == "read"
    assert records[0].sequence == "ACG"


def test_read_paired_fasta() -> None:
    path = "tests/data/simple.fasta"
    with dnaio.open(path, path) as f:
        list(f)


@pytest.mark.parametrize("interleaved", [False, True])
def test_paired_opener(fileformat, extension, interleaved) -> None:
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
        with dnaio.open(path1, path2, opener=my_opener) as f:
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


def test_detect_fastq_from_content() -> None:
    """FASTQ file that is not named .fastq"""
    with dnaio.open("tests/data/missingextension") as f:
        record = next(iter(f))
        assert record.name == "prefix:1_13_573/1"


def test_detect_compressed_fastq_from_content() -> None:
    """Compressed FASTQ file that is not named .fastq.gz"""
    with dnaio.open("tests/data/missingextension.gz") as f:
        record = next(iter(f))
    assert record.name == "prefix:1_13_573/1"


def test_detect_bam_from_content() -> None:
    with dnaio.open("tests/data/simplebamnoextension") as f:
        record = next(iter(f))
        assert record.name == "Myheader"


def test_detect_bam_from_filename() -> None:
    with dnaio.open("tests/data/simple.unaligned.bam") as f:
        record = next(iter(f))
        assert record.name == "Myheader"


def test_write(tmp_path, extension) -> None:
    out_fastq = tmp_path / ("out.fastq" + extension)
    with dnaio.open(str(out_fastq), mode="w") as f:
        f.write(dnaio.SequenceRecord("name", "ACGT", "HHHH"))
    with xopen(out_fastq) as f:
        assert f.read() == "@name\nACGT\n+\nHHHH\n"


def test_write_with_xopen(tmp_path, fileformat, extension) -> None:
    s = dnaio.SequenceRecord("name", "ACGT", "HHHH")
    out_fastq = tmp_path / ("out." + fileformat + extension)
    with xopen(out_fastq, "wb") as outer_f:
        with dnaio.open(outer_f, mode="w", fileformat=fileformat) as f:
            f.write(s)

    with xopen(out_fastq) as f:
        if fileformat == "fasta":
            assert f.read() == ">name\nACGT\n"
        else:
            assert f.read() == "@name\nACGT\n+\nHHHH\n"


def test_write_str_path(tmp_path, fileformat, extension) -> None:
    s1 = dnaio.SequenceRecord("s1", "ACGT", "HHHH")
    path = str(tmp_path / ("out." + fileformat + extension))
    with dnaio.open(path, mode="w") as f:
        f.write(s1)
    if fileformat == "fasta":
        expected = b">s1\nACGT\n"
    else:
        expected = b"@s1\nACGT\n+\nHHHH\n"
    with xopen(path, "rb") as f:
        assert f.read() == expected


def test_write_paired_same_path(tmp_path) -> None:
    path1 = tmp_path / "same.fastq"
    path2 = tmp_path / "same.fastq"
    with pytest.raises(ValueError):
        with dnaio.open(path1, path2, mode="w"):
            pass  # pragma: no cover


def test_write_paired(tmp_path, fileformat, extension) -> None:
    r1 = [
        dnaio.SequenceRecord("s1", "ACGT", "HHHH"),
        dnaio.SequenceRecord("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.SequenceRecord("t1", "TCGT", "5HHH"),
        dnaio.SequenceRecord("t2", "TGCA", "5383"),
    ]
    path1 = tmp_path / ("out.1." + fileformat + extension)
    path2 = tmp_path / ("out.2." + fileformat + extension)

    with dnaio.open(path1, path2, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
        f.write(r1[1], r2[1])
    with xopen(path1) as f:
        assert formatted_sequences(r1, fileformat) == f.read()
    with xopen(path2) as f:
        assert formatted_sequences(r2, fileformat) == f.read()


def test_write_interleaved(tmp_path, fileformat, extension) -> None:
    r1 = [
        dnaio.SequenceRecord("s1", "ACGT", "HHHH"),
        dnaio.SequenceRecord("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.SequenceRecord("t1", "TCGT", "5HHH"),
        dnaio.SequenceRecord("t2", "TGCA", "5383"),
    ]
    path = tmp_path / ("out.interleaved." + fileformat + extension)

    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
        f.write(r1[1], r2[1])
    expected = [r1[0], r2[0], r1[1], r2[1]]
    with xopen(path) as f:
        assert formatted_sequences(expected, fileformat) == f.read()


def test_append(tmp_path, fileformat, extension) -> None:
    s1 = dnaio.SequenceRecord("s1", "ACGT", "HHHH")
    s2 = dnaio.SequenceRecord("s2", "CGCA", "8383")
    path = tmp_path / ("out." + fileformat + extension)
    with dnaio.open(path, mode="w") as f:
        f.write(s1)
    with dnaio.open(path, mode="a") as f:
        f.write(s2)
    with xopen(path) as f:
        assert formatted_sequences([s1, s2], fileformat) == f.read()


def test_append_paired(tmp_path, fileformat, extension) -> None:
    r1 = [
        dnaio.SequenceRecord("s1", "ACGT", "HHHH"),
        dnaio.SequenceRecord("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.SequenceRecord("t1", "TCGT", "5HHH"),
        dnaio.SequenceRecord("t2", "TGCA", "5383"),
    ]
    path1 = tmp_path / ("out.1." + fileformat + extension)
    path2 = tmp_path / ("out.2." + fileformat + extension)

    with dnaio.open(path1, path2, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
    with dnaio.open(path1, path2, fileformat=fileformat, mode="a") as f:
        f.write(r1[1], r2[1])
    with xopen(path1) as f:
        assert formatted_sequences(r1, fileformat) == f.read()
    with xopen(path2) as f:
        assert formatted_sequences(r2, fileformat) == f.read()


def test_append_interleaved(tmp_path, fileformat, extension) -> None:
    r1 = [
        dnaio.SequenceRecord("s1", "ACGT", "HHHH"),
        dnaio.SequenceRecord("s2", "CGCA", "8383"),
    ]
    r2 = [
        dnaio.SequenceRecord("t1", "TCGT", "5HHH"),
        dnaio.SequenceRecord("t2", "TGCA", "5383"),
    ]
    path = tmp_path / ("out.interleaved." + fileformat + extension)

    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="w") as f:
        f.write(r1[0], r2[0])
    with dnaio.open(path, interleaved=True, fileformat=fileformat, mode="a") as f:
        f.write(r1[1], r2[1])
    expected = [r1[0], r2[0], r1[1], r2[1]]
    with xopen(path) as f:
        assert formatted_sequences(expected, fileformat) == f.read()


def make_random_fasta(path, n_records) -> None:
    from random import choice

    with xopen(path, "w") as f:
        for i in range(n_records):
            name = "sequence_{}".format(i)
            sequence = "".join(choice("ACGT") for _ in range(300))
            print(">", name, "\n", sequence, sep="", file=f)


def test_islice_gzip_does_not_fail(tmp_path) -> None:
    path = tmp_path / "file.fasta.gz"
    make_random_fasta(path, 100)
    f = dnaio.open(path)
    next(iter(f))
    f.close()


def test_unsupported_mode() -> None:
    with pytest.raises(ValueError) as error:
        _ = dnaio.open(os.devnull, mode="x")  # type: ignore
    error.match("Mode must be")


def test_no_file2_with_multiple_args() -> None:
    with pytest.raises(ValueError) as error:
        _ = dnaio.open(os.devnull, os.devnull, file2=os.devnull)  # type: ignore
    error.match("as positional argument")
    error.match("file2")


def test_no_multiple_files_interleaved() -> None:
    with pytest.raises(ValueError) as error:
        _ = dnaio.open(os.devnull, os.devnull, interleaved=True)  # type: ignore
    error.match("interleaved")
    error.match("one file")


@pytest.mark.parametrize(
    ["mode", "expected_class"],
    [("r", dnaio.PairedEndReader), ("w", dnaio.PairedEndWriter)],
)
def test_paired_open_with_multiple_args(
    tmp_path, fileformat, mode, expected_class
) -> None:
    path = tmp_path / "file"
    path2 = tmp_path / "file2"
    path.touch()
    path2.touch()
    with dnaio.open(path, path2, fileformat=fileformat, mode=mode) as f:
        assert isinstance(f, expected_class)


@pytest.mark.parametrize(
    ["kwargs", "expected_class"],
    [
        ({}, dnaio.multipleend.MultipleFileReader),
        ({"mode": "w"}, dnaio.multipleend.MultipleFastqWriter),
        ({"mode": "w", "fileformat": "fastq"}, dnaio.multipleend.MultipleFastqWriter),
        ({"mode": "w", "fileformat": "fasta"}, dnaio.multipleend.MultipleFastaWriter),
    ],
)
def test_multiple_open_fastq(kwargs, expected_class) -> None:
    with dnaio.open(os.devnull, os.devnull, os.devnull, **kwargs) as f:
        assert isinstance(f, expected_class)


def test_deprecated_file1_file2_keyword_arguments(tmp_path):
    path = Path("tests/data/simple.fasta")
    expected = SIMPLE_RECORDS["fasta"]
    with dnaio.open(file1=path) as f:
        records = list(f)
    assert records == expected

    with dnaio.open(path, file2=path) as f:
        records = list(f)
    assert records == list(zip(expected, expected))

    with dnaio.open(file1=path, file2=path) as f:
        records = list(f)
    assert records == list(zip(expected, expected))


def test_positional_with_file1():
    with pytest.raises(ValueError) as error:
        with dnaio.open("in.fastq", file1="in2.fastq"):
            pass  # pragma: no cover
    error.match("file1 keyword argument cannot be used together")


def test_positional_with_file1_and_file2():
    with pytest.raises(ValueError) as error:
        with dnaio.open("in.fastq", file1="in2.fastq", file2="in3.fastq"):
            pass  # pragma: no cover
    error.match("cannot be used together")
