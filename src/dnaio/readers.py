"""
Classes for reading FASTA and FASTQ files
"""
import io
from xopen import xopen
from ._sequence import Sequence
from ._util import shorten
from .exceptions import FileFormatError


class BinaryFileReader:
    """
    A mixin for readers that ensures that a file or a path can be passed in to the constructor.
    """
    _close_on_exit = False
    paired = False
    mode = 'rb'

    def __init__(self, file):
        """
        The file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).
        """
        if isinstance(file, str):
            file = xopen(file, self.mode)
            self._close_on_exit = True
        self._file = file

    def close(self):
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        if self._file is None:
            raise ValueError("I/O operation on closed BinaryFileReader")
        return self

    def __exit__(self, *args):
        self.close()


class FastaReader(BinaryFileReader):
    """
    Reader for FASTA files.
    """

    def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence):
        """
        file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).

        keep_linebreaks -- whether to keep newline characters in the sequence
        """
        super().__init__(file)
        self.sequence_class = sequence_class
        self.delivers_qualities = False
        self._delimiter = '\n' if keep_linebreaks else ''

    def __iter__(self):
        """
        Read next entry from the file (single entry at a time).
        """
        name = None
        seq = []
        f = io.TextIOWrapper(self._file)
        for i, line in enumerate(f):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == '>':
                if name is not None:
                    yield self.sequence_class(name, self._delimiter.join(seq), None)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FileFormatError("At line {0}: Expected '>' at beginning of "
                    "FASTA record, but got {1!r}.".format(i + 1, shorten(line)))

        if name is not None:
            yield self.sequence_class(name, self._delimiter.join(seq), None)
        # Prevent TextIOWrapper from closing the underlying file
        f.detach()


def _sequence_names_match(r1, r2):
    """
    Check whether the sequence records r1 and r2 have identical names, ignoring a
    suffix of '1' or '2'. Some old paired-end reads have names that end in '/1'
    and '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ)
    appends a .1 and .2 to paired-end reads if option -I is used.
    """
    name1 = r1.name.split(None, 1)[0]
    name2 = r2.name.split(None, 1)[0]
    if name1[-1:] in '12' and name2[-1:] in '12':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2


class PairedSequenceReader:
    """
    Read paired-end reads from two files.

    Wraps two BinaryFileReader instances, making sure that reads are properly
    paired.
    """
    paired = True

    def __init__(self, file1, file2, colorspace=False, fileformat=None):
        from . import open as dnaio_open  # import locally to avoid circular import
        self.reader1 = dnaio_open(file1, colorspace=colorspace, fileformat=fileformat)
        self.reader2 = dnaio_open(file2, colorspace=colorspace, fileformat=fileformat)
        self.delivers_qualities = self.reader1.delivers_qualities

    def __iter__(self):
        """
        Iterate over the paired reads. Each item is a pair of Sequence objects.
        """
        # Avoid usage of zip() below since it will consume one item too many.
        it1, it2 = iter(self.reader1), iter(self.reader2)
        while True:
            try:
                r1 = next(it1)
            except StopIteration:
                # End of file 1. Make sure that file 2 is also at end.
                try:
                    next(it2)
                    raise FileFormatError("Reads are improperly paired. There are more reads in "
                        "file 2 than in file 1.")
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FileFormatError("Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.")
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Read name '{0}' "
                    "in file 1 does not match '{1}' in file 2.".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader1.close()
        self.reader2.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceReader:
    """
    Read paired-end reads from an interleaved FASTQ file.
    """
    paired = True

    def __init__(self, file, colorspace=False, fileformat=None):
        from . import open as dnaio_open  # import locally to avoid circular import
        self.reader = dnaio_open(file, colorspace=colorspace, fileformat=fileformat)
        self.delivers_qualities = self.reader.delivers_qualities

    def __iter__(self):
        # Avoid usage of zip() below since it will consume one item too many.
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FileFormatError("Interleaved input file incomplete: Last record "
                    "{!r} has no partner.".format(r1.name))
            if not _sequence_names_match(r1, r2):
                raise FileFormatError("Reads are improperly paired. Name {0!r} "
                    "(first) does not match {1!r} (second).".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
