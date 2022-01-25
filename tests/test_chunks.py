from pytest import raises
from io import BytesIO

from dnaio._core import paired_fastq_heads
from dnaio.chunks import _fastq_head, _fasta_head, read_chunks, read_paired_chunks


def test_fasta_head():
    assert _fasta_head(b'') == 0
    assert _fasta_head(b'>1\n') == 0
    assert _fasta_head(b'>1\n3') == 0
    assert _fasta_head(b'>1\n3\n') == 0
    assert _fasta_head(b'>1\n3\n>') == 5
    assert _fasta_head(b'>1\n3\n>6') == 5
    assert _fasta_head(b'>1\n3\n>6\n') == 5
    assert _fasta_head(b'>1\n3\n>6\n8') == 5
    assert _fasta_head(b'>1\n3\n>6\n8\n') == 5
    assert _fasta_head(b'>1\n3\n>6\n8\n0') == 5
    assert _fasta_head(b'>1\n3\n>6\n8\n0\n') == 5
    assert _fasta_head(b'>1\n3\n>6\n8\n0\n>') == 12


def test_fasta_head_with_comment():
    assert _fasta_head(b'#') == 0
    assert _fasta_head(b'#\n') == 0
    assert _fasta_head(b'#\n>') == 2
    assert _fasta_head(b'#\n>3') == 2
    assert _fasta_head(b'#\n>3\n') == 2
    assert _fasta_head(b'#\n>3\n5') == 2
    assert _fasta_head(b'#\n>3\n5\n') == 2
    assert _fasta_head(b'#\n>3\n5\n>') == 7


def test_paired_fastq_heads():
    buf1 = b'first\nsecond\nthird\nfourth\nfifth'
    buf2 = b'a\nb\nc\nd\ne\nf\ng'
    assert paired_fastq_heads(buf1, buf2, len(buf1), len(buf2)) == (
        len(b'first\nsecond\nthird\nfourth\n'), len(b'a\nb\nc\nd\n'))

    assert paired_fastq_heads(b'abc', b'def', 3, 3) == (0, 0)
    assert paired_fastq_heads(b'abc\n', b'def', 4, 3) == (0, 0)
    assert paired_fastq_heads(b'abc', b'def\n', 3, 4) == (0, 0)
    assert paired_fastq_heads(b'\n\n\n\n', b'\n\n\n\n', 4, 4) == (4, 4)


def test_fastq_head():
    assert _fastq_head(b'') == 0
    assert _fastq_head(b'A\n') == 0
    assert _fastq_head(b'A\nB') == 0
    assert _fastq_head(b'A\nB\n') == 0
    assert _fastq_head(b'A\nB\nC') == 0
    assert _fastq_head(b'A\nB\nC\n') == 0
    assert _fastq_head(b'A\nB\nC\nD') == 0
    assert _fastq_head(b'A\nB\nC\nD\n') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\n') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\n') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG\n') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG\nH') == 0
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG\nH\n') == 16
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG\nH\nI') == 16
    assert _fastq_head(b'A\nB\nC\nD\nE\nF\nG\nH\nI\n') == 16


def test_read_paired_chunks():
    with open('tests/data/paired.1.fastq', 'rb') as f1:
        with open('tests/data/paired.2.fastq', 'rb') as f2:
            for c1, c2 in read_paired_chunks(f1, f2, buffer_size=128):
                print(c1, c2)


def test_read_chunks():
    for data in [b'@r1\nACG\n+\nHHH\n', b'>r1\nACGACGACG\n']:
        assert [m.tobytes() for m in read_chunks(BytesIO(data))] == [data]

        # Buffer too small
        with raises(OverflowError):
            list(read_chunks(BytesIO(data), buffer_size=4))


def test_read_chunks_empty():
    assert list(read_chunks(BytesIO(b''))) == []
