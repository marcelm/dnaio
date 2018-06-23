class FileFormatError(Exception):
    """
    The file is not formatted correctly
    """


class UnknownFileType(Exception):
    """
    The file format could not be automatically detected
    """
