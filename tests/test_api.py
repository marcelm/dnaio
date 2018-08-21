import dnaio


def test_open():
    with dnaio.open('tests/data/simple.fasta') as f:
        records = list(f)
        assert len(records) == 2


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
