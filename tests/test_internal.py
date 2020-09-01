import os
import shutil
import subprocess
import sys
from io import BytesIO
from tempfile import mkdtemp
from textwrap import dedent

from pytest import raises, mark

import dnaio
from dnaio import (
    FileFormatError, FastaFormatError, FastqFormatError,
    FastaReader, FastqReader, InterleavedSequenceReader,
    FastaWriter, FastqWriter, InterleavedSequenceWriter,
    PairedSequenceReader)
from dnaio import _record_names_match, Sequence


# files tests/data/simple.fast{q,a}
simple_fastq = [
    Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
    Sequence("second_sequence", "SEQUENCE2", "83<??:(61")
]

simple_fasta = [Sequence(x.name, x.sequence, None) for x in simple_fastq]

tiny_fastq = b'@r1\nACG\n+\nHHH\n@r2\nT\n+\n#\n'


class TestSequence:
    def test_too_many_qualities(self):
        with raises(ValueError):
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
            """# a comment
            # another one
            unexpected
            >first_sequence
            SEQUENCE1
            >second_sequence
            SEQUENCE2
            """).encode())
        with raises(FastaFormatError) as info:
            list(FastaReader(fasta))
        assert info.value.line == 2

    def test_fastareader_keeplinebreaks(self):
        with FastaReader("tests/data/simple.fasta", keep_linebreaks=True) as f:
            reads = list(f)
        assert reads[0] == simple_fasta[0]
        assert reads[1].sequence == 'SEQUEN\nCE2'

    def test_context_manager(self):
        filename = "tests/data/simple.fasta"
        with open(filename, 'rb') as f:
            assert not f.closed
            _ = list(dnaio.open(f))
            assert not f.closed
        assert f.closed

        with FastaReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            _ = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None
        # Open it a second time
        with FastaReader(filename):
            pass


class TestFastqReader:
    def test_fastqreader(self):
        with FastqReader("tests/data/simple.fastq") as f:
            reads = list(f)
        assert reads == simple_fastq

    @mark.parametrize("buffer_size", [1, 2, 3, 5, 7, 10, 20])
    def test_fastqreader_buffersize(self, buffer_size):
        with FastqReader("tests/data/simple.fastq", buffer_size=buffer_size) as f:
            reads = list(f)
        assert reads == simple_fastq

    def test_fastqreader_buffersize_too_small(self):
        with raises(ValueError):
            with FastqReader("tests/data/simple.fastq", buffer_size=0) as f:
                _ = list(f)  # pragma: no cover

    def test_fastqreader_dos(self):
        # DOS line breaks
        with open('tests/data/dos.fastq', 'rb') as f:
            assert b'\r\n' in f.read()
        with FastqReader("tests/data/dos.fastq") as f:
            dos_reads = list(f)
        with FastqReader("tests/data/small.fastq") as f:
            unix_reads = list(f)
        assert dos_reads == unix_reads

    def test_fastq_wrongformat(self):
        with raises(FastqFormatError) as info:
            with FastqReader("tests/data/withplus.fastq") as f:
                list(f)  # pragma: no cover
        assert info.value.line == 2

    def test_empty_fastq(self):
        with FastqReader(BytesIO(b'')) as fq:
            assert list(fq) == []

    @mark.parametrize('s,line', [
        (b'@', 0),
        (b'@r', 0),
        (b'@r1', 0),
        (b'@r1\n', 1),
        (b'@r1\nA', 1),
        (b'@r1\nAC', 1),
        (b'@r1\nACG', 1),
        (b'@r1\nACG\n', 2),
        (b'@r1\nACG\n+', 2),
        (b'@r1\nACG\n+\n', 3),
        (b'@r1\nACG\n+\nH', 3),
        (b'@r1\nACG\n+\nHH', 3),
        (b'@r1\nACG\n+\nHHH\n@', 4),
        (b'@r1\nACG\n+\nHHH\n@r', 4),
        (b'@r1\nACG\n+\nHHH\n@r2', 4),
        (b'@r1\nACG\n+\nHHH\n@r2\n', 5),
        (b'@r1\nACG\n+\nHHH\n@r2\nT', 5),
        (b'@r1\nACG\n+\nHHH\n@r2\nT\n', 6),
        (b'@r1\nACG\n+\nHHH\n@r2\nT\n+', 6),
        (b'@r1\nACG\n+\nHHH\n@r2\nT\n+\n', 7),
    ])
    def test_fastq_incomplete(self, s, line):
        fastq = BytesIO(s)
        with raises(FastqFormatError) as info:
            with FastqReader(fastq) as fq:
                list(fq)
        assert info.value.line == line

    def test_half_record_line_numbers(self):
        fastq = BytesIO(b'@r\nACG\n+\nHH\n')
        # Choose the buffer size such that only parts of the record fit
        # We want to ensure that the line number is reset properly
        # after the record has been half-parsed
        buffer_size = len('@r\nACG\n+\n')
        with raises(FastqFormatError) as info:
            with FastqReader(fastq, buffer_size=buffer_size) as fq:
                list(fq)  # pragma: no cover
        assert 'Length of sequence and qualities differ' in info.value.message
        assert info.value.line == 3

    @mark.parametrize('s,line', [
        (b'@r1\nACG\n+\nH#HH\n@r2\nT\n+\nH\n', 3),
        (b'@r1\nACG\n+\n#H\n@r2\nT\n+\nH\n', 3),
        (b'@r1\nACG\n+\nHHH\n@r2\nT\n+\nHH\n', 7),
        (b'@r1\nACG\n+\nHHH\n@r2\nT\n+\n\n', 7),
    ])
    def test_differing_lengths(self, s, line):
        fastq = BytesIO(s)
        with raises(FastqFormatError) as info:
            with FastqReader(fastq) as fq:
                list(fq)
        assert info.value.line == line

    def test_missing_final_newline(self):
        # Files with a missing final newline are currently allowed
        fastq = BytesIO(b'@r1\nA\n+\nH')
        with dnaio.open(fastq) as f:
            records = list(f)
        assert records == [Sequence('r1', 'A', 'H')]

    def test_not_opened_as_binary(self):
        filename = 'tests/data/simple.fastq'
        with open(filename, 'rt') as f:
            with raises(ValueError):
                list(dnaio.open(f))

    def test_context_manager(self):
        filename = "tests/data/simple.fastq"
        with open(filename, 'rb') as f:
            assert not f.closed
            _ = list(dnaio.open(f))
            assert not f.closed
        assert f.closed

        with FastqReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            _ = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None

    def test_two_header_detection(self):
        fastq = BytesIO(b'@r1\nACG\n+r1\nHHH\n@r2\nT\n+r2\n#\n')
        with FastqReader(fastq) as fq:
            assert fq.two_headers
            list(fq)

        fastq = BytesIO(b'@r1\nACG\n+\nHHH\n@r2\nT\n+r2\n#\n')
        with FastqReader(fastq) as fq:
            assert not fq.two_headers
            list(fq)

    def test_second_header_not_equal(self):
        fastq = BytesIO(b'@r1\nACG\n+xy\n')
        with raises(FastqFormatError) as info:
            with FastqReader(fastq) as fq:
                list(fq)  # pragma: no cover
        assert "Sequence descriptions don't match" in info.value.message


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

    def test_autodetect_fasta_format(self, tmpdir):
        path = str(tmpdir.join('tmp.fasta'))
        with dnaio.open(path, mode='w') as f:
            assert isinstance(f, FastaWriter)
            for seq in simple_fastq:
                f.write(seq)
        with dnaio.open(path) as f:
            records = list(f)
        assert records == simple_fasta

    def test_write_qualities_to_fasta(self):
        path = os.path.join(self._tmpdir, 'tmp.fasta')
        with dnaio.open(path, mode='w', qualities=True) as f:
            assert isinstance(f, FastaWriter)
            for seq in simple_fastq:
                f.write(seq)
        with dnaio.open(path) as f:
            assert list(f) == simple_fasta

    def test_autodetect_fastq_format(self):
        path = os.path.join(self._tmpdir, 'tmp.fastq')
        with dnaio.open(path, mode='w') as f:
            assert isinstance(f, FastqWriter)
            for seq in simple_fastq:
                f.write(seq)
        with dnaio.open(path) as f:
            assert list(f) == simple_fastq

    def test_autodetect_fastq_weird_name(self):
        path = os.path.join(self._tmpdir, 'tmp.fastq.gz')
        with dnaio.open(path, mode='w') as f:
            assert isinstance(f, FastqWriter)
            for seq in simple_fastq:
                f.write(seq)
        weird_path = os.path.join(self._tmpdir, 'tmp.weird.gz')
        os.rename(path, weird_path)
        with dnaio.open(weird_path) as f:
            assert list(f) == simple_fastq

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
        with InterleavedSequenceReader("tests/data/interleaved.fastq") as isr:
            reads = list(isr)

        assert reads == expected
        with dnaio.open("tests/data/interleaved.fastq", interleaved=True) as f:
            reads = list(f)
        assert reads == expected

    def test_missing_partner(self):
        s = BytesIO(b'@r1\nACG\n+\nHHH\n')
        with raises(FileFormatError) as info:
            with InterleavedSequenceReader(s) as isr:
                list(isr)
        assert "Interleaved input file incomplete" in info.value.message

    def test_incorrectly_paired(self):
        s = BytesIO(b'@r1/1\nACG\n+\nHHH\n@wrong_name\nTTT\n+\nHHH\n')
        with raises(FileFormatError) as info:
            with InterleavedSequenceReader(s) as isr:
                list(isr)
        assert "Reads are improperly paired" in info.value.message


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
        assert bio.getvalue() == b'>name\n\n', '{!r}'.format(bio.getvalue())


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
        assert bio.getvalue() == (
            b'@A/1 comment\nTTA\n+\n##H\n'
            b'@A/2 comment\nGCT\n+\nHH#\n'
            b'@B/1\nCC\n+\nHH\n'
            b'@B/2\nTG\n+\n#H\n'
        )


class TestPairedSequenceReader:
    def test_read(self):
        s1 = BytesIO(b'@r1\nACG\n+\nHHH\n')
        s2 = BytesIO(b'@r2\nGTT\n+\n858\n')
        with PairedSequenceReader(s1, s2) as psr:
            assert [
                (Sequence("r1", "ACG", "HHH"), Sequence("r2", "GTT", "858")),
            ] == list(psr)

    def test_record_names_match(self):
        match = _record_names_match
        assert match('abc', 'abc')
        assert match('abc/1', 'abc/2')
        assert match('abc.1', 'abc.2')
        assert match('abc1', 'abc2')
        assert not match('abc', 'xyz')

    def test_missing_partner1(self):
        s1 = BytesIO(b'')
        s2 = BytesIO(b'@r1\nACG\n+\nHHH\n')

        with raises(FileFormatError) as info:
            with PairedSequenceReader(s1, s2) as psr:
                list(psr)
        assert "There are more reads in file 2 than in file 1" in info.value.message

    def test_missing_partner2(self):
        s1 = BytesIO(b'@r1\nACG\n+\nHHH\n')
        s2 = BytesIO(b'')

        with raises(FileFormatError) as info:
            with PairedSequenceReader(s1, s2) as psr:
                list(psr)
        assert "There are more reads in file 1 than in file 2" in info.value.message

    def test_incorrectly_paired(self):
        s1 = BytesIO(b'@r1/1\nACG\n+\nHHH\n')
        s2 = BytesIO(b'@wrong_name\nTTT\n+\nHHH\n')
        with raises(FileFormatError) as info:
            with PairedSequenceReader(s1, s2) as psr:
                list(psr)
        assert "Reads are improperly paired" in info.value.message


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
    py = subprocess.Popen(
        [sys.executable, 'tests/read_from_stdin.py'], stdin=cat.stdout, stdout=subprocess.PIPE)
    cat.stdout.close()

    # Check that the read_from_stdin.py script prints the correct number of records
    assert str(expected) == py.communicate()[0].decode().strip()
