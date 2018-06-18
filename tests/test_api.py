import dnaio


def test_open():
    with dnaio.open('tests/data/simple.fasta') as f:
        pass
