from setuptools import setup, Extension
from Cython.Build import cythonize
import setuptools_scm  # noqa  Ensure it’s installed

setup(
    ext_modules=cythonize(
        [
            Extension("dnaio._core", sources=["src/dnaio/_core.pyx"]),
        ]
    ) + [
        Extension("dnaio._sequence_bytes", sources=["src/dnaio/_sequence_bytes.c"])
    ],
)
