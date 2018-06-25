from pytest import raises

from dnaio.exceptions import FileFormatError
from dnaio.colorspace import ColorspaceSequence


def test_too_many_qualities_colorspace():
    with raises(FileFormatError):
        ColorspaceSequence(name="name", sequence="T0123", qualities="#####")


def test_invalid_primer():
    with raises(FileFormatError):
        ColorspaceSequence(name="name", sequence="K0123", qualities="####")
