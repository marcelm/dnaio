=========
Changelog
=========

v1.2.0 (2023-12-11)
-------------------

* :pr:`124`: Added support for chunking FASTA reads to ``read_paired_chunks``.
  Previously, only FASTQ was supported.

v1.1.0 (2023-11-20)
-------------------

* :pr:`116`: Added experimental support for reading unaligned BAM files
  (single-end only at the moment). This uses a custom, highly efficient
  BAM parser, making dnaio faster than htslib in this particular case.

v1.0.1 (2023-10-06)
-------------------

* :pr:`120`: Improved type annotations.
* Dropped support for Python 3.7
* Added support for Python 3.12

v1.0.0 (2023-09-06)
-------------------
* :pr:`110`: Added ``id`` and ``comment`` properties to ``SequenceRecord``.

v0.10.0 (2022-12-05)
--------------------

* :pr:`99`: SequenceRecord initialization is now faster, which also provides
  a speed boost to FASTQ iteration. ``SequenceRecord.__new__`` cannot be used
  anymore to initialize `SequenceRecord` objects.
* :pr:`96`: ``open_threads`` and ``compression_level`` are now added
  to `~dnaio.open` as arguments. By default dnaio now uses compression level
  1 and does not utilize external programs to speed up gzip (de)compression.
* :pr:`87`: `~dnaio.open` can now open more than two files.
  The ``file1`` and ``file2`` arguments are now deprecated.

v0.9.1 (2022-08-01)
-------------------

* :pr:`85`: macOS wheels are now also built as part of the release procedure.
* :pr:`81`: API documentation improvements and minor code refactors for
  readability.

v0.9.0 (2022-05-17)
-------------------

* :pr:`79`: Added a `~dnaio.records_are_mates` function to be used for checking whether
  three or more records are mates of each other (by checking the ID).
* :pr:`74`, :pr:`68`: Made FASTQ parsing faster by implementing the check for
  ASCII using SSE vector instructions.
* :pr:`72`: Added a `tutorial <https://dnaio.readthedocs.io/en/latest/tutorial.html>`_.

v0.8.0 (2022-03-26)
-------------------

* Preliminary documentation is available at
  <https://dnaio.readthedocs.io/>.
* :pr:`53`: Renamed ``Sequence`` to `~dnaio.SequenceRecord`.
  The previous name is still available as an alias
  so that existing code will continue to work.
* When reading a FASTQ file, there is now a check that ensures that
  all characters are ASCII.
* Function ``record_names_match`` is deprecated, use `~dnaio.SequenceRecord.is_mate` instead.
* Added `~dnaio.SequenceRecord.reverse_complement`.
* Dropped Python 3.6 support as it is end-of-life.

v0.7.1 (2022-01-26)
-------------------

* :pr:`34`: Fix parsing of FASTA files that just contain a comment and no reads

v0.7.0 (2022-01-17)
-------------------

* @rhpvorderman contributed many performance improvements in :pr:`15`,
  :pr:`17`, :pr:`18`, :pr:`20`, :pr:`21`, :pr:`22`, :pr:`23`. Reading
  and writing FASTQ files and reading of paired-end FASTQ files was
  sped up significantly. For example, reading uncompressed FASTQ is
  50% faster (!) than before.
* :pr:`28`: Windows support added


v0.6.0 (2021-09-28)
-------------------

* :pr:`12`: Improve FASTQ writing speed twofold (thanks to @rhpvorderman)


v0.5.2 (2021-09-07)
-------------------

* :issue:`7`: Ignore a trailing "3" in the read id
