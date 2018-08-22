class FileFormatError(Exception):
    """
    The file is not formatted correctly
    """


class UnknownFileFormat(Exception):
    """
    The file format could not be automatically detected
    """
