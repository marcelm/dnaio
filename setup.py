import os
import platform
import sys

from setuptools import setup, Extension
import setuptools_scm  # noqa  Ensure it’s installed

if platform.machine() == "AMD64":
    # Macro is defined by default for clang and GCC on relevant targets, but
    # not by MSVC.
    DEFINE_MACROS = [("__SSE2__", 1)]
else:
    DEFINE_MACROS = []


def extra_compile_args():
    if sys.platform.startswith("linux") and platform.machine() == "x86_64":
        # DNAIO_CPU_BASIC=1 pip install --no-binary dnaio dnaio
        # Will work for linux users that do not have SSSE3 support.
        if not os.getenv("DNAIO_CPU_BASIC"):
            return ["-mssse3"]
    # Do not bother with Windows and MacOS as it is not given that they have
    # a compiler installed. Simply use compatible wheels instead so no users
    # run into trouble.
    return None


setup(
    ext_modules=[
        Extension(
            "dnaio._core",
            sources=["src/dnaio/_core.pyx"],
            define_macros=DEFINE_MACROS,
            extra_compile_args=extra_compile_args(),
        ),
    ],
)
