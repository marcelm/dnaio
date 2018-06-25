from .exceptions import FileFormatError
from ._util import shorten
from ._sequence import Sequence
from ._core import FastqReader
from .readers import FastaReader
from .writers import FastaWriter, FastqWriter


class ColorspaceSequence(Sequence):
    def __init__(self, name, sequence, qualities, primer=None):
        # In colorspace, the first character is the last nucleotide of the primer base
        # and the second character encodes the transition from the primer base to the
        # first real base of the read.
        if primer is None:
            self.primer = sequence[0:1]
            sequence = sequence[1:]
        else:
            self.primer = primer
        if qualities is not None and len(sequence) != len(qualities):
            rname = shorten(name)
            raise FileFormatError("In read named {0!r}: length of colorspace quality "
                "sequence ({1}) and length of read ({2}) do not match (primer "
                "is: {3!r})".format(rname, len(qualities), len(sequence), self.primer))
        super().__init__(name, sequence, qualities)
        if self.primer not in ('A', 'C', 'G', 'T'):
            raise FileFormatError("Primer base is {0!r} in read {1!r}, but it "
                "should be one of A, C, G, T.".format(
                    self.primer, shorten(name)))

    def __repr__(self):
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={0!r}'.format(shorten(self.qualities))
        return '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'.format(
            shorten(self.name), self.primer, shorten(self.sequence), qstr)

    def __getitem__(self, key):
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.qualities is not None else None,
            self.primer)

    def __reduce__(self):
        return (ColorspaceSequence, (self.name, self.sequence, self.qualities, self.primer,
            self.second_header, self.match))


def sra_colorspace_sequence(name, sequence, qualities, second_header):
    """Factory for an SRA colorspace sequence (which has one quality value too many)"""
    return ColorspaceSequence(name, sequence, qualities[1:])


class FastaQualReader:
    """
    Reader for reads that are stored in .(CS)FASTA and .QUAL files.
    """
    delivers_qualities = True
    paired = False

    def __init__(self, fastafile, qualfile, sequence_class=Sequence):
        """
        fastafile and qualfile are filenames or file-like objects.
        If a filename is used, then .gz files are recognized.

        The objects returned when iteritng over this file are instances of the
        given sequence_class.
        """
        self.fastareader = FastaReader(fastafile)
        self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
        self.sequence_class = sequence_class

    def __iter__(self):
        """
        Yield Sequence objects.
        """
        # conversion dictionary: maps strings to the appropriate ASCII-encoded character
        conv = dict()
        for i in range(-5, 256 - 33):
            conv[str(i)] = chr(i + 33)
        for fastaread, qualread in zip(self.fastareader, self.qualreader):
            if fastaread.name != qualread.name:
                raise FileFormatError("The read names in the FASTA and QUAL file "
                    "do not match ({0!r} != {1!r})".format(fastaread.name, qualread.name))
            try:
                qualities = ''.join([conv[value] for value in qualread.sequence.split()])
            except KeyError as e:
                raise FileFormatError("Within read named {0!r}: Found invalid quality "
                    "value {1}".format(fastaread.name, e))
            assert fastaread.name == qualread.name
            yield self.sequence_class(fastaread.name, fastaread.sequence, qualities)

    def close(self):
        self.fastareader.close()
        self.qualreader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ColorspaceFastaReader(FastaReader):
    def __init__(self, file, keep_linebreaks=False):
        super().__init__(file, keep_linebreaks, sequence_class=ColorspaceSequence)


class ColorspaceFastqReader(FastqReader):
    def __init__(self, file):
        super().__init__(file, sequence_class=ColorspaceSequence)


class SRAColorspaceFastqReader(FastqReader):
    def __init__(self, file):
        super().__init__(file, sequence_class=sra_colorspace_sequence)


class ColorspaceFastaQualReader(FastaQualReader):
    def __init__(self, fastafile, qualfile):
        super().__init__(fastafile, qualfile, sequence_class=ColorspaceSequence)


class ColorspaceFastaWriter(FastaWriter):
    def write(self, record):
        name = record.name
        sequence = record.primer + record.sequence
        super().write(name, sequence)


class ColorspaceFastqWriter(FastqWriter):
    def write(self, record):
        name = record.name
        sequence = record.primer + record.sequence
        qualities = record.qualities
        super().writeseq(name, sequence, qualities)
