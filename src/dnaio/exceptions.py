from typing import Optional


class FileFormatError(Exception):
    """
    The file is not formatted correctly
    """
    format = 'sequence'  # Something generic that works for both FASTA and FASTQ

    def __init__(self, msg: str, line: Optional[int]):
        super().__init__(msg, line)
        self.message = msg
        self.line = line  # starts at 0!

    def __str__(self):
        line = 'unknown line' if self.line is None else 'line {}'.format(self.line + 1)
        return 'Error in {} file at {}: {}'.format(self.format, line, self.message)


class FastqFormatError(FileFormatError):
    format = 'FASTQ'


class FastaFormatError(FileFormatError):
    format = 'FASTA'


class UnknownFileFormat(Exception):
    """
    The file format could not be automatically detected
    """
