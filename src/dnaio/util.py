from xopen import xopen


class FormatError(Exception):
    """
    Raised when an input file is malformatted.
    """


def _shorten(s, n=100):
    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[:n-3] + '...'
    return s


class BinaryFileReader:
    """Read possibly compressed files containing sequences"""
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
