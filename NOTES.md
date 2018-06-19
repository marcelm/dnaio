
- compressed
- paired-end
- interleaved
- chunked
- FASTA, FASTQ
- BAM?
-


import dnaio

with dnaio.open('input.fastq.gz') as f:
    for record in f:
        print(record....)

with dnaio.open('input.1.fastq.gz', 'input.2.fastq.gz') as f:
    for record in f:
        print(record....)



Use cases

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
- Should SequenceRecord be immutable?


TODO

- Sequence.name should be Sequence.description or so (reserve .name for the part
  before the first space)

- FormatErrors should have a line_number attribute
Documentation

- Line endings
- second header
- multi-line FASTQ


FASTQ chunks

- need an index attribute
- need a line_number attribute
