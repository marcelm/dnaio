import textwrap

from pytest import raises
from io import BytesIO

import dnaio
from dnaio import UnknownFileFormat, FileFormatError
from dnaio._core import paired_fastq_heads
from dnaio.chunks import (
    _fastq_head,
    _fasta_head,
    read_chunks,
    read_paired_chunks,
    _paired_fasta_heads,
)


def test_fasta_head():
    assert _fasta_head(b"") == 0
    assert _fasta_head(b">1\n") == 0
    assert _fasta_head(b">1\n3") == 0
    assert _fasta_head(b">1\n3\n") == 0
    assert _fasta_head(b">1\n3\n>") == 5
    assert _fasta_head(b">1\n3\n>6") == 5
    assert _fasta_head(b">1\n3\n>6\n") == 5
    assert _fasta_head(b">1\n3\n>6\n8") == 5
    assert _fasta_head(b">1\n3\n>6\n8\n") == 5
    assert _fasta_head(b">1\n3\n>6\n8\n0") == 5
    assert _fasta_head(b">1\n3\n>6\n8\n0\n") == 5
    assert _fasta_head(b">1\n3\n>6\n8\n0\n>") == 12


def test_fasta_head_with_comment():
    assert _fasta_head(b"#") == 0
    assert _fasta_head(b"#\n") == 0
    assert _fasta_head(b"#\n>") == 2
    assert _fasta_head(b"#\n>3") == 2
    assert _fasta_head(b"#\n>3\n") == 2
    assert _fasta_head(b"#\n>3\n5") == 2
    assert _fasta_head(b"#\n>3\n5\n") == 2
    assert _fasta_head(b"#\n>3\n5\n>") == 7


def test_paired_fasta_heads():
    def pheads(buf1, buf2):
        return _paired_fasta_heads(buf1, buf2, len(buf1), len(buf2))

    assert pheads(b"", b"") == (0, 0)
    assert pheads(b">r", b">r") == (0, 0)
    assert pheads(b">r\nA\n>s", b">r") == (0, 0)
    assert pheads(b">r\nA\n>s", b">r\nCT\n>s") == (5, 6)
    assert pheads(b">r\nA\n>s\nG\n>t\n", b">r\nCT\n>s") == (5, 6)

    buf1 = (
        textwrap.dedent(
            """
        >1
        a
        b
        >2
        c
        >3
        uv
        """
        )
        .strip()
        .encode()
    )
    buf2 = (
        textwrap.dedent(
            """
        >1
        def
        >2
        gh
        i
        >3
        """
        )
        .strip()
        .encode()
    )

    assert pheads(buf1, buf2) == (
        len(b">1\na\nb\n>2\nc\n"),
        len(b">1\ndef\n>2\ngh\ni\n"),
    )


def test_paired_fastq_heads():
    buf1 = b"first\nsecond\nthird\nfourth\nfifth"
    buf2 = b"a\nb\nc\nd\ne\nf\ng"
    assert paired_fastq_heads(buf1, buf2, len(buf1), len(buf2)) == (
        len(b"first\nsecond\nthird\nfourth\n"),
        len(b"a\nb\nc\nd\n"),
    )

    assert paired_fastq_heads(b"abc", b"def", 3, 3) == (0, 0)
    assert paired_fastq_heads(b"abc\n", b"def", 4, 3) == (0, 0)
    assert paired_fastq_heads(b"abc", b"def\n", 3, 4) == (0, 0)
    assert paired_fastq_heads(b"\n\n\n\n", b"\n\n\n\n", 4, 4) == (4, 4)


def test_fastq_head():
    assert _fastq_head(b"") == 0
    assert _fastq_head(b"A\n") == 0
    assert _fastq_head(b"A\nB") == 0
    assert _fastq_head(b"A\nB\n") == 0
    assert _fastq_head(b"A\nB\nC") == 0
    assert _fastq_head(b"A\nB\nC\n") == 0
    assert _fastq_head(b"A\nB\nC\nD") == 0
    assert _fastq_head(b"A\nB\nC\nD\n") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\n") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\n") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG\n") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG\nH") == 0
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG\nH\n") == 16
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG\nH\nI") == 16
    assert _fastq_head(b"A\nB\nC\nD\nE\nF\nG\nH\nI\n") == 16


def test_read_paired_chunks_fastq():
    with open("tests/data/paired.1.fastq", "rb") as f1:
        with open("tests/data/paired.2.fastq", "rb") as f2:
            for c1, c2 in read_paired_chunks(f1, f2, buffer_size=128):
                print(c1, c2)


def test_paired_chunks_fasta(tmp_path):
    for i in (1, 2):
        with dnaio.open(f"tests/data/paired.{i}.fastq") as infile:
            with dnaio.open(tmp_path / f"{i}.fasta", mode="w") as outfile:
                for record in infile:
                    record.qualities = None
                    outfile.write(record)

    with open(tmp_path / "1.fasta", "rb") as r1:
        with open(tmp_path / "2.fasta", "rb") as r2:
            for c1, c2 in read_paired_chunks(r1, r2, buffer_size=128):
                print(c1.tobytes(), c2.tobytes())


def test_paired_chunks_different_number_of_records():
    record = b"@r\nAA\n+\n##\n"
    buf1 = record
    buf2 = record * 3
    it = read_paired_chunks(BytesIO(buf1), BytesIO(buf2), 16)
    assert next(it) == (record, record)
    with raises(FileFormatError) as error:
        next(it)
    error.match("more data found in the other file")


def test_read_chunks():
    for data in [b"@r1\nACG\n+\nHHH\n", b">r1\nACGACGACG\n"]:
        assert [m.tobytes() for m in read_chunks(BytesIO(data))] == [data]

        # Buffer too small
        with raises(OverflowError):
            list(read_chunks(BytesIO(data), buffer_size=4))


def test_read_chunks_empty():
    assert list(read_chunks(BytesIO(b""))) == []


def test_invalid_file_format():
    with raises(UnknownFileFormat):
        list(read_chunks(BytesIO(b"invalid format")))
