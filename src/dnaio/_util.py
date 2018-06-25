def shorten(s, n=100):
    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[:n-3] + '...'
    return s
