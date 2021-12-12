import os.path
from setuptools import setup, Extension, find_packages
from setuptools.command.sdist import sdist
from setuptools.command.build_ext import build_ext


def no_cythonize(extensions, **_ignore):
    """Change .pyx to .c or .cpp (copied from Cython documentation)"""
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources


class BuildExt(build_ext):
    def run(self):
        # If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
        # file retrieved from PyPI that already includes the pre-cythonized
        # extension modules, and then we do not need to run cythonize().
        if os.path.exists('PKG-INFO'):
            no_cythonize(self.extensions)
        else:
            # Otherwise, this is a 'developer copy' of the code, and then the
            # only sensible thing is to require Cython to be installed.
            from Cython.Build import cythonize
            self.extensions = cythonize(self.extensions)
            # Setuptools build_meta requires a _needs_stub attribute.
            # This is required for `pip install .`
            for i in range(len(self.extensions)):
                setattr(self.extensions[i], "_needs_stub", False)
        super().run()


class SDist(sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(self.distribution.ext_modules)
        super().run()


with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='dnaio',
    setup_requires=['setuptools_scm'],  # Support pip versions that don't know about pyproject.toml
    use_scm_version={'write_to': 'src/dnaio/_version.py'},
    author='Marcel Martin',
    author_email='marcel.martin@scilifelab.se',
    url='https://github.com/marcelm/dnaio/',
    description='Read and write FASTA and FASTQ files efficiently',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={"dnaio": ["py.typed", "*.pyi"]},
    extras_require={
        'dev': ['Cython', 'pytest'],
    },
    ext_modules=[
        Extension('dnaio._core', sources=['src/dnaio/_core.pyx']),
    ],
    cmdclass={'build_ext': BuildExt, 'sdist': SDist},
    install_requires=['xopen>=0.8.2'],
    python_requires='>=3.6',
    classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Cython",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
