=========
Changelog
=========

development version
-------------------

* Renamed ``Sequence`` to ``SequenceRecord`` and ``BytesSequence`` to ``BytesSequenceRecord``.
  The previous names are still available as aliases so existing code will continue to work.
* Function ``record_names_match`` is deprecated, use `SequenceRecord.is_mate` instead.

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
