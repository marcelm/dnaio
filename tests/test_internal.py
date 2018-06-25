import os
import shutil
import subprocess
import sys
from io import BytesIO
from tempfile import mkdtemp
from textwrap import dedent

import dnaio
from dnaio import (
    two_fastq_heads, _fastq_head,
    read_paired_chunks, read_chunks_from_file, Sequence, FileFormatError,
    FastaReader, FastqReader, FastaQualReader, InterleavedSequenceReader,
    FastaWriter, FastqWriter, InterleavedSequenceWriter)

from dnaio.readers import _sequence_names_match
from pytest import raises, mark

# files tests/data/simple.fast{q,a}
simple_fastq = [
    Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
    Sequence("second_sequence", "SEQUENCE2", "83<??:(61")
]

simple_fasta = [Sequence(x.name, x.sequence, None) for x in simple_fastq]

tiny_fastq = b'@r1\nACG\n+\nHHH\n@r2\nT\n+\n#\n'


class TestSequence:
    def test_too_many_qualities(self):
        with raises(FileFormatError):
            Sequence(name="name", sequence="ACGT", qualities="#####")


class TestFastaReader:
    def test_file(self):
        with FastaReader("tests/data/simple.fasta") as f:
            reads = list(f)
        assert reads == simple_fasta

    def test_bytesio(self):
        fasta = BytesIO(b">first_sequence\nSEQUENCE1\n>second_sequence\nSEQUENCE2\n")
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_with_comments(self):
        fasta = BytesIO(dedent(
            """
            # a comment
            # another one
            >first_sequence
            SEQUENCE1
            >second_sequence
            SEQUENCE2
            """).encode())
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_wrong_format(self):
        fasta = BytesIO(dedent(
            """
            # a comment
            # another one
            unexpected
            >first_sequence
            SEQUENCE1
            >second_sequence
            SEQUENCE2
            """).encode())
        with raises(FileFormatError):
            reads = list(FastaReader(fasta))

    def test_fastareader_keeplinebreaks(self):
        with FastaReader("tests/data/simple.fasta", keep_linebreaks=True) as f:
            reads = list(f)
        assert reads[0] == simple_fasta[0]
        assert reads[1].sequence == 'SEQUEN\nCE2'

    def test_context_manager(self):
        filename = "tests/data/simple.fasta"
        with open(filename, 'rb') as f:
            assert not f.closed
            reads = list(dnaio.open(f))
            assert not f.closed
        assert f.closed

        with FastaReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            reads = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None
        # Open it a second time
        with FastaReader(filename) as sr:
            pass


class TestFastqReader:
    def test_fastqreader(self):
        with FastqReader("tests/data/simple.fastq") as f:
            reads = list(f)
        assert reads == simple_fastq

    def test_fastqreader_dos(self):
        with open('tests/data/dos.fastq', 'rb') as f:
            assert b'\r\n' in f.read()
        with FastqReader("tests/data/dos.fastq") as f:
            dos_reads = list(f)
        with FastqReader("tests/data/small.fastq") as f:
            unix_reads = list(f)
        assert dos_reads == unix_reads

    def test_fastq_wrongformat(self):
        with FastqReader("tests/data/withplus.fastq") as f:
            with raises(FileFormatError):
                reads = list(f)

    def test_empty_fastq(self):
        with FastqReader(BytesIO(b'')) as fq:
            assert list(fq) == []

    @mark.parametrize('s', [
        b'@',
        b'@r',
        b'@r1',
        b'@r1\n',
        b'@r1\nA',
        b'@r1\nAC',
        b'@r1\nACG',
        b'@r1\nACG\n',
        b'@r1\nACG\n+',
        b'@r1\nACG\n+\n',
        b'@r1\nACG\n+\nH',
        b'@r1\nACG\n+\nHH',
        b'@r1\nACG\n+\nHHH\n@',
        b'@r1\nACG\n+\nHHH\n@r',
        b'@r1\nACG\n+\nHHH\n@r2',
        b'@r1\nACG\n+\nHHH\n@r2\n',
        b'@r1\nACG\n+\nHHH\n@r2\nT',
        b'@r1\nACG\n+\nHHH\n@r2\nT\n',
        b'@r1\nACG\n+\nHHH\n@r2\nT\n+',
        b'@r1\nACG\n+\nHHH\n@r2\nT\n+\n',
    ])
    def test_fastq_incomplete(self, s):
        fastq = BytesIO(s)
        with FastqReader(fastq) as fq:
            with raises(FileFormatError):
                list(fq)

    def test_missing_final_newline(self):
        fastq = BytesIO(b'@r1\nA\n+\nH')
        with FastqReader(fastq) as fq:
            assert len(list(fq)) == 1

    def test_not_opened_as_binary(self):
        filename = 'tests/data/simple.fastq'
        with open(filename, 'rt') as f:
            with raises(ValueError):
                list(dnaio.open(f))

    def test_context_manager(self):
        filename = "tests/data/simple.fastq"
        with open(filename, 'rb') as f:
            assert not f.closed
            reads = list(dnaio.open(f))
            assert not f.closed
        assert f.closed

        with FastqReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            reads = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None


class TestFastaQualReader:
    def test_mismatching_read_names(self):
        fasta = BytesIO(b">name\nACG")
        qual = BytesIO(b">nome\n3 5 7")
        with raises(FileFormatError):
            list(FastaQualReader(fasta, qual))

    def test_invalid_quality_value(self):
        fasta = BytesIO(b">name\nACG")
        qual = BytesIO(b">name\n3 xx 7")
        with raises(FileFormatError):
            list(FastaQualReader(fasta, qual))


class TestOpen:
    def setup(self):
        self._tmpdir = mkdtemp()

    def teardown(self):
        shutil.rmtree(self._tmpdir)

    def test_sequence_reader(self):
        # test the autodetection
        with dnaio.open("tests/data/simple.fastq") as f:
            reads = list(f)
        assert reads == simple_fastq

        with dnaio.open("tests/data/simple.fasta") as f:
            reads = list(f)
        assert reads == simple_fasta

        with open("tests/data/simple.fastq", 'rb') as f:
            reads = list(dnaio.open(f))
        assert reads == simple_fastq

        # make the name attribute unavailable
        f = BytesIO(open("tests/data/simple.fastq", 'rb').read())
        reads = list(dnaio.open(f))
        assert reads == simple_fastq

        f = BytesIO(open("tests/data/simple.fasta", 'rb').read())
        reads = list(dnaio.open(f))
        assert reads == simple_fasta

    def test_autodetect_fasta_format(self):
        path = os.path.join(self._tmpdir, 'tmp.fasta')
        with dnaio.open(path, mode='w') as f:
            assert isinstance(f, FastaWriter)
            for seq in simple_fastq:
                f.write(seq)
        assert list(dnaio.open(path)) == simple_fasta

    def test_write_qualities_to_fasta(self):
        path = os.path.join(self._tmpdir, 'tmp.fasta')
        with dnaio.open(path, mode='w', qualities=True) as f:
            assert isinstance(f, FastaWriter)
            for seq in simple_fastq:
                f.write(seq)
        assert list(dnaio.open(path)) == simple_fasta

    def test_autodetect_fastq_format(self):
        path = os.path.join(self._tmpdir, 'tmp.fastq')
        with dnaio.open(path, mode='w') as f:
            assert isinstance(f, FastqWriter)
            for seq in simple_fastq:
                f.write(seq)
        assert list(dnaio.open(path)) == simple_fastq

    def test_fastq_qualities_missing(self):
        path = os.path.join(self._tmpdir, 'tmp.fastq')
        with raises(ValueError):
            dnaio.open(path, mode='w', qualities=False)


class TestInterleavedReader:
    def test(self):
        expected = [
            (
                Sequence('read1/1 some text', 'TTATTTGTCTCCAGC', '##HHHHHHHHHHHHH'),
                Sequence('read1/2 other text', 'GCTGGAGACAAATAA', 'HHHHHHHHHHHHHHH')
            ),
            (
                Sequence('read3/1', 'CCAACTTGATATTAATAACA', 'HHHHHHHHHHHHHHHHHHHH'),
                Sequence('read3/2', 'TGTTATTAATATCAAGTTGG', '#HHHHHHHHHHHHHHHHHHH')
            ),
        ]
        reads = list(InterleavedSequenceReader("tests/data/interleaved.fastq"))
        for (r1, r2), (e1, e2) in zip(reads, expected):
            print(r1, r2, e1, e2)

        assert reads == expected
        with dnaio.open("tests/data/interleaved.fastq", interleaved=True) as f:
            reads = list(f)
        assert reads == expected

    def test_missing_partner(self):
        s = BytesIO(b'@r1\nACG\n+\nHHH')
        with raises(FileFormatError):
            list(InterleavedSequenceReader(s))

    def test_incorrectly_paired(self):
        s = BytesIO(b'@r1/1\nACG\n+\nHHH\n@wrong_name\nTTT\n+\nHHH')
        with raises(FileFormatError):
            list(InterleavedSequenceReader(s))


class TestFastaWriter:
    def setup(self):
        self._tmpdir = mkdtemp()
        self.path = os.path.join(self._tmpdir, 'tmp.fasta')

    def teardown(self):
        shutil.rmtree(self._tmpdir)

    def test(self):
        with FastaWriter(self.path) as fw:
            fw.write("name", "CCATA")
            fw.write("name2", "HELLO")
        assert fw._file.closed
        with open(self.path) as t:
            assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

    def test_linelength(self):
        with FastaWriter(self.path, line_length=3) as fw:
            fw.write("r1", "ACG")
            fw.write("r2", "CCAT")
            fw.write("r3", "TACCAG")
        assert fw._file.closed
        with open(self.path) as t:
            d = t.read()
            assert d == '>r1\nACG\n>r2\nCCA\nT\n>r3\nTAC\nCAG\n'

    def test_write_sequence_object(self):
        with FastaWriter(self.path) as fw:
            fw.write(Sequence("name", "CCATA"))
            fw.write(Sequence("name2", "HELLO"))
        assert fw._file.closed
        with open(self.path) as t:
            assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

    def test_write_to_file_like_object(self):
        bio = BytesIO()
        with FastaWriter(bio) as fw:
            fw.write(Sequence("name", "CCATA"))
            fw.write(Sequence("name2", "HELLO"))
        assert bio.getvalue() == b'>name\nCCATA\n>name2\nHELLO\n'
        assert not bio.closed
        assert not fw._file.closed

    def test_write_zero_length_sequence(self):
        bio = BytesIO()
        with FastaWriter(bio) as fw:
            fw.write(Sequence("name", ""))
        assert bio.getvalue() == b'>name\n\n', '{0!r}'.format(bio.getvalue())


class TestFastqWriter:
    def setup(self):
        self._tmpdir = mkdtemp()
        self.path = os.path.join(self._tmpdir, 'tmp.fastq')

    def teardown(self):
        shutil.rmtree(self._tmpdir)

    def test(self):
        with FastqWriter(self.path) as fq:
            fq.writeseq("name", "CCATA", "!#!#!")
            fq.writeseq("name2", "HELLO", "&&&!&&")
        assert fq._file.closed
        with open(self.path) as t:
            assert t.read() == '@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'

    def test_twoheaders(self):
        with FastqWriter(self.path, two_headers=True) as fq:
            fq.write(Sequence("name", "CCATA", "!#!#!"))
            fq.write(Sequence("name2", "HELLO", "&&&!&"))
        assert fq._file.closed
        with open(self.path) as t:
            assert t.read() == '@name\nCCATA\n+name\n!#!#!\n@name2\nHELLO\n+name2\n&&&!&\n'

    def test_write_to_file_like_object(self):
        bio = BytesIO()
        with FastqWriter(bio) as fq:
            fq.writeseq("name", "CCATA", "!#!#!")
            fq.writeseq("name2", "HELLO", "&&&!&&")
        assert bio.getvalue() == b'@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'


class TestInterleavedWriter:
    def test(self):
        reads = [
            (
                Sequence('A/1 comment', 'TTA', '##H'),
                Sequence('A/2 comment', 'GCT', 'HH#')
            ),
            (
                Sequence('B/1', 'CC', 'HH'),
                Sequence('B/2', 'TG', '#H')
            ),
        ]
        bio = BytesIO()
        with InterleavedSequenceWriter(bio) as writer:
            for read1, read2 in reads:
                writer.write(read1, read2)
        assert bio.getvalue() == b'@A/1 comment\nTTA\n+\n##H\n@A/2 comment\nGCT\n+\nHH#\n@B/1\nCC\n+\nHH\n@B/2\nTG\n+\n#H\n'


class TestPairedSequenceReader:
    def test_sequence_names_match(self):
        def match(name1, name2):
            seq1 = Sequence(name1, 'ACGT')
            seq2 = Sequence(name2, 'AACC')
            return _sequence_names_match(seq1, seq2)

        assert match('abc', 'abc')
        assert match('abc/1', 'abc/2')
        assert match('abc.1', 'abc.2')
        assert match('abc1', 'abc2')
        assert not match('abc', 'xyz')


def test_two_fastq_heads():
    buf1 = b'first\nsecond\nthird\nfourth\nfifth'
    buf2 = b'a\nb\nc\nd\ne\nf\ng'
    assert two_fastq_heads(buf1, buf2, len(buf1), len(buf2)) == (
        len(b'first\nsecond\nthird\nfourth\n'), len(b'a\nb\nc\nd\n'))

    assert two_fastq_heads(b'abc', b'def', 3, 3) == (0, 0)
    assert two_fastq_heads(b'abc\n', b'def', 4, 3) == (0, 0)
    assert two_fastq_heads(b'abc', b'def\n', 3, 4) == (0, 0)
    assert two_fastq_heads(b'\n\n\n\n', b'\n\n\n\n', 4, 4) == (4, 4)


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


def test_read_chunks_from_file():
    for data in [b'@r1\nACG\n+\nHHH\n', b'>r1\nACGACGACG\n']:
        assert [m.tobytes() for m in read_chunks_from_file(BytesIO(data))] == [data]

        # Buffer too small
        with raises(OverflowError):
            list(read_chunks_from_file(BytesIO(data), buffer_size=4))


@mark.parametrize('path', [
    'tests/data/simple.fastq',
    'tests/data/dos.fastq',
    'tests/data/simple.fasta',
    'tests/data/with_comment.fasta',
])
def test_read_stdin(path):
    # Get number of records in the input file
    with dnaio.open(path) as f:
        expected = len(list(f))

    # Use 'cat' to simulate that no file name is available for stdin of the subprocess
    cat = subprocess.Popen(['cat', path], stdout=subprocess.PIPE)
    py = subprocess.Popen([sys.executable, 'tests/read_from_stdin.py'], stdin=cat.stdout,
        stdout=subprocess.PIPE)
    cat.stdout.close()

    # Check that the read_from_stdin.py script prints the correct number of records
    assert str(expected) == py.communicate()[0].decode().strip()
