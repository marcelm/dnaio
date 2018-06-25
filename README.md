# dnaio parses FASTQ and FASTA

`dnaio` is a Python 3 library for fast parsing of FASTQ and also FASTA files. The code was previously part of the
[cutadapt](https://cutadapt.readthedocs.io/) tool and has been improved since it has been split out.


## Example usage

The main interface is the `dnaio.open` function ::

    import dnaio

    with dnaio.open('reads.fastq.gz') as f:
        bp = 0
        for record in f:
            bp += len(record)
    print(f'The input file contains {bp/1E6:.1f} Mbp')


## Features and supported file types

- FASTQ input and output
- FASTA input and output
- Compressed input and output (`.gz`, `.bz2` and `.xz`)
- Paired-end data in two files
- Interleaved paired-end data in a single file
- Files with DOS/Windows linebreaks can be read
- FASTQ files with a second header line (after the `+`) are supported


# Limitations

- Multi-line FASTQ files are not supported. You shouldn’t use them anyway.
- FASTQ parsing is the focus of this library. The FASTA parser is not as optimized.


# Links

* [Source code](https://github.com/marcelm/dnaio/)
* [Report an issue](https://github.com/marcelm/dnaio/issues)
* [Project page on PyPI](https://pypi.python.org/pypi/dnaio/)
