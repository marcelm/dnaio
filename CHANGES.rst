=========
Changelog
=========

v0.7.1 (2022-01-26)
-------------------

* PR #34: Fix parsing of FASTA files that just contain a comment and no reads

v0.7.0 (2022-01-17)
-------------------

* @rhpvorderman contributed many performance improvements in PR #15, #17, #18, #20, #21, #22, #23. Reading and writing FASTQ files and reading of paired-end FASTQ files was sped up significantly. For example, reading uncompressed FASTQ is 50% faster (!) than before.
* PR #28: Windows support added


v0.6.0 (2021-09-28)
-------------------

* PR #12: Improve FASTQ writing speed twofold (thanks to @rhpvorderman)


v0.5.2 (2021-09-07)
-------------------

* Issue #7: Ignore a trailing "3" in the read id
