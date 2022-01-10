from setuptools import setup, Extension
from Cython.Build import cythonize
import setuptools_scm  # noqa  Ensure it’s installed

setup(
    ext_modules=cythonize(
        [
            Extension("dnaio._sequence", sources=["src/dnaio/_sequence.pyx"]),
            Extension("dnaio._core", sources=["src/dnaio/_core.pyx"]),
        ]
    ),
)
