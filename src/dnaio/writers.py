import io
from xopen import xopen


class FileWriter:
    file_mode = 'wt'

    def __init__(self, file, _close_file=None):
        assert self.file_mode in ('wt', 'wb')
        if isinstance(file, str):
            self._file = xopen(file, self.file_mode)
            self._close_on_exit = True
        else:
            if self.file_mode == 'wb':
                self._file = self._force_binary_stream(file)
            else:
                self._file = self._force_text_stream(file)
            self._close_on_exit = bool(_close_file)

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

    def __init__(self, file, line_length=None, _close_file=None):
        """
        If line_length is not None, the lines will
        be wrapped after line_length characters.
        """
        super().__init__(file, _close_file=_close_file)
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
            if not sequence:
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
    file_mode = 'wb'

    def __init__(self, file, two_headers=False, _close_file=None):
        super().__init__(file, _close_file=_close_file)
        self._two_headers = two_headers

    def write(self, record):
        """
        Write a Sequence record to the FASTQ file.

        The record object must have attributes .name, .sequence and .qualities.
        """
        name2 = record.name if self._two_headers else ''
        s = ('@' + record.name + '\n' + record.sequence + '\n+'
             + name2 + '\n' + record.qualities + '\n')
        self._file.write(s.encode('ascii'))

    def writeseq(self, name, sequence, qualities):
        self._file.write("@{0:s}\n{1:s}\n+\n{2:s}\n".format(
            name, sequence, qualities).encode('ascii'))
