.. image:: https://github.com/marcelm/dnaio/workflows/CI/badge.svg
    :alt: GitHub Actions badge

.. image:: https://img.shields.io/pypi/v/dnaio.svg?branch=main
    :target: https://pypi.python.org/pypi/dnaio
    :alt: PyPI badge

.. image:: https://codecov.io/gh/marcelm/dnaio/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/marcelm/dnaio
    :alt: Codecov badge

=====================================
dnaio processes FASTQ and FASTA files
=====================================

``dnaio`` is a Python 3.8+ library for very efficient parsing and writing of FASTQ and also FASTA files.
The code was previously part of the
`Cutadapt <https://cutadapt.readthedocs.io/>`_ tool and has been improved significantly since it has been split out.

Example usage
=============

The main interface is the `dnaio.open <https://dnaio.readthedocs.io/en/latest/api.html>`_ function::

    import dnaio

    with dnaio.open("reads.fastq.gz") as f:
        bp = 0
        for record in f:
            bp += len(record)
    print(f"The input file contains {bp/1E6:.1f} Mbp")

For more, see the `tutorial <https://dnaio.readthedocs.io/en/latest/tutorial.html>`_ and
`API documentation <https://dnaio.readthedocs.io/en/latest/api.html>`_.

Installation
============

Using pip:: 

    pip install dnaio zstandard

``zstandard`` can be omitted if support for Zstandard (``.zst``) files is not required.

Features and supported file types
=================================

- FASTQ input and output
- FASTA input and output
- BAM input
- Compressed input and output (``.gz``, ``.bz2``, ``.xz`` and ``.zst`` are detected automatically)
- Paired-end data in two files
- Interleaved paired-end data in a single file
- Files with DOS/Windows linebreaks can be read
- FASTQ files with a second header line (after the ``+``) are supported

Limitations
===========

- Multi-line FASTQ files are not supported
- FASTQ and BAM parsing is the focus of this library. The FASTA parser is not as optimized

Links
=====

* `Documentation <https://dnaio.readthedocs.io/>`_
* `Source code <https://github.com/marcelm/dnaio/>`_
* `Report an issue <https://github.com/marcelm/dnaio/issues>`_
* `Project page on PyPI <https://pypi.python.org/pypi/dnaio/>`_
