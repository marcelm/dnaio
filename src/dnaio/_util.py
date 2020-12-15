import pathlib


def _is_path(obj: object) -> bool:
    """
    Return whether the given object looks like a path (str, pathlib.Path or pathlib2.Path)
    """
    # pytest uses pathlib2.Path objects on Python 3.5 for its tmp_path fixture.
    # On Python 3.6+, this function can be replaced with isinstance(obj, os.PathLike)
    import sys
    if "pathlib2" in sys.modules:
        import pathlib2  # type: ignore
        path_classes = [str, pathlib.Path, pathlib2.Path]
    else:
        path_classes = [str, pathlib.Path]
    return isinstance(obj, tuple(path_classes))


def shorten(s: str, n: int = 100) -> str:

    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[:n-3] + '...'
    return s
