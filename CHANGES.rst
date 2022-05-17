=========
Changelog
=========

v0.9.0 (2022-05-17)
-------------------

* :pr:`79`: Added a `records_are_mates` function to be used for checking whether
  three or more records are mates of each other (by checking the ID).
* :pr:`74`, :pr:`68`: Made FASTQ parsing faster by implementing the check for
  ASCII using SSE vector instructions.
* :pr:`72`: Added a `tutorial <https://dnaio.readthedocs.io/en/latest/tutorial.html>`_.

v0.8.0 (2022-03-26)
-------------------

* Preliminary documentation is available at
  <https://dnaio.readthedocs.io/>.
* :pr:`53`: Renamed ``Sequence`` to `SequenceRecord`.
  The previous name is still available as an alias
  so that existing code will continue to work.
* When reading a FASTQ file, there is now a check that ensures that
  all characters are ASCII.
* Function ``record_names_match`` is deprecated, use `SequenceRecord.is_mate` instead.
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
