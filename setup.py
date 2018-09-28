import sys
import os.path
from setuptools import setup, Extension
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext
import versioneer

if sys.version_info[:2] < (3, 4):
    sys.stdout.write('Python 3.4 or later is required\n')
    sys.exit(1)


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


extensions = [
    Extension('dnaio._core', sources=['src/dnaio/_core.pyx']),
]

cmdclass = versioneer.get_cmdclass()
versioneer_build_ext = cmdclass.get('build_ext', _build_ext)
versioneer_sdist = cmdclass.get('sdist', _sdist)


class build_ext(versioneer_build_ext):
    def run(self):
        # If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
        # file retrieved from PyPI that already includes the pre-cythonized
        # extension modules, and then we do not need to run cythonize().
        if os.path.exists('PKG-INFO'):
            no_cythonize(extensions)
        else:
            # Otherwise, this is a 'developer copy' of the code, and then the
            # only sensible thing is to require Cython to be installed.
            from Cython.Build import cythonize
            self.extensions = cythonize(self.extensions)
        versioneer_build_ext.run(self)


class sdist(versioneer_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(extensions)
        versioneer_sdist.run(self)


cmdclass['build_ext'] = build_ext
cmdclass['sdist'] = sdist

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='dnaio',
    version=versioneer.get_version(),
    author='Marcel Martin',
    author_email='marcel.martin@scilifelab.se',
    url='https://github.com/marcelm/dnaio/',
    description='Read FASTA and FASTQ files efficiently',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    packages=['dnaio'],
    package_dir={'': 'src'},
    extras_require={
        'dev': ['Cython', 'pytest'],
    },
    ext_modules=extensions,
    cmdclass=cmdclass,
    install_requires=['xopen'],
    classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Cython",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
