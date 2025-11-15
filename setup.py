import platform

from setuptools import setup, Extension
import setuptools_scm  # noqa  Ensure it’s installed

if platform.machine() == "AMD64":
    # Macro is defined by default for clang and GCC on relevant targets, but
    # not by MSVC.
    DEFINE_MACROS = [("__SSE2__", 1)]
else:
    DEFINE_MACROS = []

DEFINE_MACROS += [
    ("CYTHON_LIMITED_API", 1),
    ("Py_Limited_API", 0x030C0000),
]

setup(
    ext_modules=[
        Extension(
            "dnaio._core",
            sources=["src/dnaio/_core.pyx"],
            define_macros=DEFINE_MACROS,
            py_limited_api=True,
        ),
    ],
)
