import pathlib


def _is_path(obj: object) -> bool:
    """
    Return whether the given object looks like a path (str or pathlib.Path)
    """
    return isinstance(obj, (str, pathlib.Path))


def shorten(s: str, n: int = 100) -> str:

    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[: n - 3] + "..."
    return s
