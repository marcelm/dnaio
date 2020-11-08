from dnaio._util import shorten


def test_shorten():
    assert shorten(None) is None
    assert shorten("hello too long", 5) == "he..."
    assert shorten("hello not too long") == "hello not too long"
