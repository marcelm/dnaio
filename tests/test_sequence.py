import sys

import pytest

from dnaio import SequenceRecord, BytesSequenceRecord


class TestSequence:
    def test_too_many_qualities(self):
        with pytest.raises(ValueError):
            SequenceRecord(name="name", sequence="ACGT", qualities="#####")

    def test_fastq_bytes(self):
        assert SequenceRecord("name", "ACGT", "====").fastq_bytes() == \
            b"@name\nACGT\n+\n====\n"

    def test_fastq_bytes_two_headers(self):
        assert SequenceRecord("name", "ACGT", "====").fastq_bytes_two_headers() == \
            b"@name\nACGT\n+name\n====\n"

    def test_is_mate_succes(self):
        assert SequenceRecord("name1", "A", "=").is_mate(SequenceRecord("name2", "GC", "FF"))

    def test_init_name_bad(self):
        with pytest.raises(ValueError) as error:
            SequenceRecord("nąme1", "A", "=")
        error.match("ASCII")

    def test_init_name_none(self):
        with pytest.raises(TypeError) as error:
            SequenceRecord(None, "A", "=")
        error.match("str")

    def test_init_sequence_bad(self):
        with pytest.raises(ValueError) as error:
            SequenceRecord("name1", "Ä", "=")
        error.match("ASCII")

    def test_init_sequence_none(self):
        with pytest.raises(TypeError) as error:
            SequenceRecord("name1", None, "=")
        error.match("str")

    def test_init_qualities_bad(self):
        with pytest.raises(ValueError) as error:
            SequenceRecord("name1", "A", "ä")
        error.match("ASCII")

    def test_init_qualities_none(self):
        seq = SequenceRecord("name1", "A", None)
        assert seq.qualities is None

    def test_set_name_bad(self):
        seq = SequenceRecord("name1", "A", "=")
        with pytest.raises(ValueError) as error:
            seq.name = "näme1"
        error.match("ASCII")

    def test_set_name_none(self):
        seq = SequenceRecord("name1", "A", "=")
        with pytest.raises(TypeError) as error:
            seq.name = None
        error.match("str")

    def test_set_sequence_bad(self):
        seq = SequenceRecord("name1", "A", "=")
        with pytest.raises(ValueError) as error:
            seq.sequence = "Ä"
        error.match("ASCII")

    def test_set_sequence_none(self):
        seq = SequenceRecord("name1", "A", "=")
        with pytest.raises(TypeError) as error:
            seq.sequence = None
        error.match("str")

    def test_set_qualities_bad(self):
        seq = SequenceRecord("name1", "A", "=")
        with pytest.raises(ValueError) as error:
            seq.qualities = "Ä"
        error.match("ASCII")

    def test_set_qualities_none(self):
        seq = SequenceRecord("name1", "A", "=")
        seq.qualities = None
        assert seq.qualities is None


class TestBytesSequence:
    def test_too_many_qualities(self):
        with pytest.raises(ValueError):
            BytesSequenceRecord(name=b"name", sequence=b"ACGT", qualities=b"#####")

    def test_fastq_bytes(self):
        assert BytesSequenceRecord(b"name", b"ACGT", b"====").fastq_bytes() == \
            b"@name\nACGT\n+\n====\n"

    def test_fastq_bytes_two_headers(self):
        seq = BytesSequenceRecord(b"", b"", b"")
        # Below creates an invalid sequence, but this is done to see if the
        # underlying function properly takes into account lengths of the
        # attributes.
        seq.name = b"name"
        seq.sequence = b"ACGTA"
        seq.qualities = b"=="
        assert seq.fastq_bytes_two_headers() == b"@name\nACGTA\n+name\n==\n"

    def test_reference_counts(self):
        # Make sure BytesSequence is properly implemented so there are no
        # reference leaks.
        name = b"name"
        sequence = b"ACGT"
        qualities = b"===="
        name_ref = sys.getrefcount(name)
        seq_ref = sys.getrefcount(sequence)
        qual_ref = sys.getrefcount(qualities)
        seqbytes = BytesSequenceRecord(name, sequence, qualities)
        assert sys.getrefcount(name) == name_ref + 1
        assert sys.getrefcount(sequence) == seq_ref + 1
        assert sys.getrefcount(qualities) == qual_ref + 1
        del seqbytes
        assert sys.getrefcount(name) == name_ref
        assert sys.getrefcount(sequence) == seq_ref
        assert sys.getrefcount(qualities) == qual_ref

    def test_is_mate_succes(self):
        assert BytesSequenceRecord(b"name1", b"A", b"=").is_mate(
            BytesSequenceRecord(b"name2", b"GC", b"FF"))

    def test_init_name_none(self):
        with pytest.raises(TypeError) as error:
            BytesSequenceRecord(None, b"A", b"=")
        error.match("bytes")

    def test_init_sequence_none(self):
        with pytest.raises(TypeError) as error:
            BytesSequenceRecord(b"name1", None, b"=")
        error.match("bytes")

    def test_init_qualities_none(self):
        seq = BytesSequenceRecord(b"name1", b"A", None)
        assert seq.qualities is None

    def test_init_qualities_wrong_tpye(self):
        with pytest.raises(TypeError) as error:
            BytesSequenceRecord(b"name1", b"A", "=")
        error.match("bytes")

    def test_set_name_none(self):
        seq = BytesSequenceRecord(b"name1", b"A", b"=")
        with pytest.raises(TypeError) as error:
            seq.name = None
        error.match("bytes")

    def test_set_sequence_none(self):
        seq = BytesSequenceRecord(b"name1", b"A", b"=")
        with pytest.raises(TypeError) as error:
            seq.sequence = None
        error.match("bytes")

    def test_set_qualities_none(self):
        seq = BytesSequenceRecord(b"name1", b"A", b"=")
        seq.qualities = None
        assert seq.qualities is None
