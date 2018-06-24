
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
- optimize writing
- FileFormatErrors should have a line_number attribute
Documentation

- test FASTQ reader with smaller buffer sizes
- Line endings
- second header
- multi-line FASTQ
- accept FASTQ files that are missing the final linebreak

FASTQ chunks

- need an index attribute
- need a line_number attribute


# API

## Advertised

- dnaio.open
- Sequence(Record)
- possibly SequencePair/PairedSequence?


## Reader

- FastqReader
- FastaReader
- PairedSequenceReader -> rename to PairedFastqReader?
- InterleavedSequenceReader -> rename to InterleavedFastqReader


## Writing

class FastqWriter
class FastaWriter
class PairedSequenceWriter
class InterleavedSequenceWriter


## Chunking

def find_fasta_record_end(buf, end):
def find_fastq_record_end(buf, end=None):
def read_chunks_from_file(f, buffer_size=4*1024**2):
def read_paired_chunks(f, f2, buffer_size=4*1024**2):
head
fastq_head
two_fastq_heads


## Colorspace

class ColorspaceSequence(Sequence):
def sra_colorspace_sequence(name, sequence, qualities, second_header):
class ColorspaceFastaReader(FastaReader):
class ColorspaceFastqReader(FastqReader):
class SRAColorspaceFastqReader(FastqReader):
class ColorspaceFastaQualReader(FastaQualReader):
class ColorspaceFastaWriter(FastaWriter):
- FastaQualReader
class ColorspaceFastqWriter
