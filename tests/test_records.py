import pytest

from dnaio import SequenceRecord


class TestSequenceRecord:
    def test_too_many_qualities(self):
        with pytest.raises(ValueError):
            SequenceRecord(name="name", sequence="ACGT", qualities="#####")

    def test_fastq_bytes(self):
        assert (
            SequenceRecord("name", "ACGT", "====").fastq_bytes()
            == b"@name\nACGT\n+\n====\n"
        )

    def test_fastq_bytes_two_headers(self):
        assert (
            SequenceRecord("name", "ACGT", "====").fastq_bytes(two_headers=True)
            == b"@name\nACGT\n+name\n====\n"
        )

    def test_is_mate_succes(self):
        assert SequenceRecord("name1", "A", "=").is_mate(
            SequenceRecord("name2", "GC", "FF")
        )

    def test_reverse_complement(self):
        assert SequenceRecord(
            "name1",
            "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
            "/AAAA/6E/EEEEEEEEEEEE/EEEEA///E/",
        ).reverse_complement() == SequenceRecord(
            "name1",
            "nvhdbmrswykaacgtNVHDBMRSWYKAACGT",
            "/E///AEEEE/EEEEEEEEEEEE/E6/AAAA/",
        )

    def test_reverse_complement_none_qualities(self):
        assert SequenceRecord(
            "name1", "GATTACA", None
        ).reverse_complement() == SequenceRecord("name1", "TGTAATC", None)

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

    def test_set_id(self):
        seq = SequenceRecord("name", "A", "=")
        with pytest.raises(AttributeError):
            seq.id = "Obi-Wan"

    def test_set_comment(self):
        seq = SequenceRecord("name", "A", "=")
        with pytest.raises(AttributeError):
            seq.comment = "Hello there!"

    @pytest.mark.parametrize(
        ["record", "expected"],
        [
            (SequenceRecord("name", "A", "="), None),
            (SequenceRecord("name ", "A", "="), None),
            (SequenceRecord("name  ", "A", "="), None),
            (SequenceRecord("name", "A", "="), None),
            (SequenceRecord("AotC I hate sand!", "A", "="), "I hate sand!"),
            (
                SequenceRecord("Givemesome                       space", "A", "="),
                "space",
            ),
        ],
    )
    def test_get_comment(self, record, expected):
        assert record.comment == expected

    @pytest.mark.parametrize(
        ["record", "expected"],
        [
            (SequenceRecord("name", "A", "="), "name"),
            (SequenceRecord("name ", "A", "="), "name"),
            (SequenceRecord("name  ", "A", "="), "name"),
            (SequenceRecord("name", "A", "="), "name"),
            (SequenceRecord("AotC I hate sand!", "A", "="), "AotC"),
        ],
    )
    def test_get_id(self, record, expected):
        assert record.id == expected

    def test_reset_id_and_comment_on_name_update(self):
        record = SequenceRecord("Obi-Wan: don't try it!", "", "")
        assert record.id == "Obi-Wan:"
        assert record.comment == "don't try it!"
        record.name = "Anakin: you underestimate my power!"
        assert record.id == "Anakin:"
        assert record.comment == "you underestimate my power!"


def test_legacy_sequence():
    from dnaio import Sequence

    s = Sequence("name", "ACGT", "####")
    assert isinstance(s, SequenceRecord)
