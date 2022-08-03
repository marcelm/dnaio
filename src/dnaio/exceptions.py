from typing import Optional


class UnknownFileFormat(Exception):
    """
    The file format could not be automatically detected
    """


class FileFormatError(Exception):
    """
    The file is not formatted correctly

    Attributes:
        line: If available, the number of the line at which the error occurred
          or None if not.
          The first line has index 0.
    """

    format = "sequence"  # Something generic that works for both FASTA and FASTQ

    def __init__(self, msg: str, line: Optional[int]):
        super().__init__(msg, line)
        self.message = msg
        self.line = line  # starts at 0!

    def __str__(self):
        line = "unknown line" if self.line is None else f"line {self.line + 1}"
        return f"Error in {self.format} file at {line}: {self.message}"


class FastqFormatError(FileFormatError):
    """
    The FASTQ file is not formatted correctly
    """

    format = "FASTQ"


class FastaFormatError(FileFormatError):
    """
    The FASTA file is not formatted correctly
    """

    format = "FASTA"
