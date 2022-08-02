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


class IsNotAscii(ValueError):
    """
    Some part of a FASTA or FASTQ file is not ASCII encoded

    Attributes:
        field: Description of the field ("name", "sequence", "qualities" or similar)
            in which non-ASCII characters were found
        value: Unicode string that was intended to be assigned to the field
    """

    def __init__(self, field: str, value: str):
        self.value = value
        self.field = field

    def __str__(self):
        detail = ""
        try:
            self.value.encode("ascii")
        except UnicodeEncodeError as e:
            detail = (
                f", but found '{self.value[e.start:e.end]}' at position {e.start + 1}"
            )
        return f"'{self.field}' in sequence file must be ASCII encoded{detail}"
