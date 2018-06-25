import io
from xopen import xopen


class FileWriter:
    file_mode = 'wt'

    def __init__(self, file):
        assert self.file_mode in ('wt', 'wb')
        if isinstance(file, str):
            self._file = xopen(file, self.file_mode)
            self._close_on_exit = True
        else:
            if self.file_mode == 'wb':
                self._file = self._force_binary_stream(file)
            else:
                self._file = self._force_text_stream(file)
            self._close_on_exit = False

    @staticmethod
    def _force_binary_stream(file):
        if hasattr(file, 'readinto'):
            return file
        elif isinstance(file, io.StringIO):
            return io.BytesIO(file.getvalue().encode('ascii'))
        else:
            return file.buffer

    @staticmethod
    def _force_text_stream(file):
        if hasattr(file, 'readinto'):
            import codecs
            return codecs.getwriter('ascii')(file)
        else:
            return file

    def close(self):
        if self._close_on_exit:
            self._file.close()

    def __enter__(self):
        if self._file.closed:
            raise ValueError("I/O operation on closed file")
        return self

    def __exit__(self, *args):
        self.close()


class FastaWriter(FileWriter):
    """
    Write FASTA-formatted sequences to a file.
    """

    def __init__(self, file, line_length=None):
        """
        If line_length is not None, the lines will
        be wrapped after line_length characters.
        """
        super().__init__(file)
        self.line_length = line_length if line_length != 0 else None

    def write(self, name_or_record, sequence=None):
        """Write an entry to the the FASTA file.

        If only one parameter (name_or_record) is given, it must have
        attributes .name and .sequence, which are then used.
        Otherwise, the first parameter must be the name and the second
        the sequence.

        The effect is that you can write this:
        writer.write("name", "ACCAT")
        or
        writer.write(Sequence("name", "ACCAT"))
        """
        if sequence is None:
            name = name_or_record.name
            sequence = name_or_record.sequence
        else:
            name = name_or_record

        if self.line_length is not None:
            print('>{0}'.format(name), file=self._file)
            for i in range(0, len(sequence), self.line_length):
                print(sequence[i:i + self.line_length], file=self._file)
            if len(sequence) == 0:
                print(file=self._file)
        else:
            print('>{0}'.format(name), sequence, file=self._file, sep='\n')


class FastqWriter(FileWriter):
    """
    Write sequences with qualities in FASTQ format.

    FASTQ files are formatted like this:
    @read name
    SEQUENCE
    +
    QUALITIS
    """

    def __init__(self, file, two_headers=False):
        super().__init__(file)
        self._two_headers = two_headers

    def write(self, record):
        """
        Write a Sequence record to the the FASTQ file.

        The record must have attributes .name, .sequence and .qualities.
        """
        name2 = record.name if self._two_headers else ''
        s = ('@' + record.name + '\n' + record.sequence + '\n+'
             + name2 + '\n' + record.qualities + '\n')
        self._file.write(s)

    def writeseq(self, name, sequence, qualities):
        print("@{0:s}\n{1:s}\n+\n{2:s}".format(
            name, sequence, qualities), file=self._file)


class PairRecordWriter:
    """Public interface to paired-record files"""

    def write(self, read1, read2):
        raise NotImplementedError()

    def close(self):
        raise NotImplementedError()

    def __enter__(self):
        # TODO do not allow this twice
        return self

    def __exit__(self, *args):
        self.close()


class PairedSequenceWriter(PairRecordWriter):
    def __init__(self, file1, file2, colorspace=False, fileformat='fastq', qualities=None):
        from . import open as dnaio_open  # import locally to avoid circular import

        self._writer1 = dnaio_open(file1, colorspace=colorspace, fileformat=fileformat, mode='w',
                             qualities=qualities)
        self._writer2 = dnaio_open(file2, colorspace=colorspace, fileformat=fileformat, mode='w',
                             qualities=qualities)

    def write(self, read1, read2):
        self._writer1.write(read1)
        self._writer2.write(read2)

    def close(self):
        self._writer1.close()
        self._writer2.close()


class InterleavedSequenceWriter(PairRecordWriter):
    """
    Write paired-end reads to an interleaved FASTA or FASTQ file
    """

    def __init__(self, file, colorspace=False, fileformat='fastq', qualities=None):
        from . import open as dnaio_open  # import locally to avoid circular import

        self._writer = dnaio_open(
            file, colorspace=colorspace, fileformat=fileformat, mode='w', qualities=qualities)

    def write(self, read1, read2):
        self._writer.write(read1)
        self._writer.write(read2)

    def close(self):
        self._writer.close()


