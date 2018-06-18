import dnaio


def test_open():
    with dnaio.open('in.fasta') as f:
        pass


"""
TODO

Test the following use cases

- open FASTQ from path
- open FASTA from path
- open compressed FASTA or FASTQ (.gz, .bz2, .xz)
- open paired-end data
- open interleaved data
- open file-like object (such as sys.stdin)
- use custom sequence record class
- autodetect file format from contents
- write FASTQ
- write FASTA
- read FASTQ/FASTA chunks (multiple records)

Issues

- Binary vs text
- FASTQ second header (after +)
- 
"""
