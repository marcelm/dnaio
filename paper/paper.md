---
title: 'dnaio: Really fast FASTQ parsing'
tags:
  - Python
  - bioinformatics
  - FASTQ
authors:
  - name: Marcel Martin
    orcid: https://orcid.org/0000-0002-0680-200X
    equal-contrib: true
    affiliation: 1
  - name: Ruben H P Vorderman
    orcid: https://orcid.org/0000-0002-8813-1528
    equal-contrib: true
    affiliation: 2
affiliations:
 - name:  Department of Biochemistry and Biophysics, National Bioinformatics Infrastructure Sweden, Science for Life Laboratory, Stockholm University, Solna, Sweden
   index: 1
 - name:  Sequencing Analysis Support Core, Department of Biomedical Data Sciences, Leiden University Medical Center, 2333 ZA, Leiden, The Netherlands
   index: 2
date: 99 Month 2022
bibliography: paper.bib

---

Your paper should include:

...

- Mention (if applicable) a representative set of past or ongoing
research projects using the software and recent scholarly publications enabled by it.


# Summary

<!--
A summary describing the high-level functionality and purpose of the
software for a diverse, non-specialist audience.
-->
Low-cost, high-throughput DNA sequencing has become the workhorse of life
science research. Compressed FASTQ files [@Cock2010Sanger] are the nearly
universally used file format for representing the nucleotide sequences
obtained from a DNA sequencing machine. Parsing FASTQ files is thus an
essential task performed in many bioinformatics pipelines. We describe
here our software `dnaio`, a highly optimized Python library
for reading and writing FASTQ files that is faster than the fastest C/C++
implementation known to us.

Given the size
of typical FASTQ files (~200 GiB after decompression for a human genome
sequenced at

dnaio is a Python library for reading and writing of FASTQ
files. FASTQ files are used to , a file format used in bioinformatics to represent the
sequences obtained from a sequencing machine ("reads").

dnaio is highly optimized
with support for
FASTQ and FASTA files are simple text-based file formats used for
storing DNA sequences. DNA sequencers typically output the
obtained sequences ("reads") in gzip-compressed FASTQ format, which
contain also nucleotide-level sequencing quality information in
addition to the nucleotide sequences themselves [@Cock2010Sanger].
Bioinformatics tools that work with FASTQ is used for

...

# Statement of need

Although alternatives have been proposed (uncompressed BAM, etc. TODO),
gzip-compressed FASTQ is still the dominant format for representing sequencing
reads with base-level quality values. It is used as the base-called output
format from sequencing machines, as a long-term storage format (TODO SRA),
but also as an intermediate format between stages in a bioinformatics pipeline.
In the latter case, disk I/O and compression/decompression overhead may be
avoided by working in a streaming manner, that is, by connecting the
FASTQ-manipulating processes via Unix pipes, but the parsing/serializing
overhead remains.

dnaio was conceived because tools such as Cutadapt `[@Martin2011Cutadapt]`
or fastq-filter (TODO) can be so fast in individual sequence manipulations
that FASTQ parsing and writing overhead dominates the runtime.
dnaio is now used in Cutadapt for FASTQ and FASTA input and output.
Because Cutadapt is downloaded more than 1000 times each day (counting
downloads from PyPI and Anaconda.org), we can assume that also dnaio, as one of
its dependencies, works well in practice on many different types of datasets.



# The library

dnaio is a Python library,


<!--
A Statement of need section that clearly illustrates the research
purpose of the software and places it in the context of related work.
-->

# To Do

* test coverage
* serves as an example that Python libraries do not need to be slow
* benchmarks and comparison to other tools
* features compared to other tools (feature matrix)


FASTQ parsing consists of a few simple steps.
- Decompressing the FASTQ file.
- Finding newlines to delimit the strings in the record.
- Check if header lines start with the correct symbol.
- Validating whether the strings contain characters in the acceptable range.
- Create Python objects in memory to represent name, sequence and qualities.

The parsing code is written in the C-code generator Cython `[@Behnel2011Cython]`
to make full use of Python's C-API. Dnaio accelerates gzip compression and
decompression by utilizing ISA-L `[@isa-l]`. Newlines are found using
the C standard library function `memchr`. Because it is a standard, it has been
well optimized. The x864-64 glibc implementation utilizes SIMD instructions for
example. All characters in FASTQ records should be in the ASCII range. This can be
checked by reducing all characters using bitwise
OR and checking if the most significant bit of the result is set. This
is sped up by using SIMD instructions. This is a trick utilized by simd-utf8
`[CITATION NEEDED]`. This algorithm is run on the entire input buffer rather
than on individual strings to further enhance performance. More extensive checking
such as checking for IUPAC only characters in the sequence string is not
performed as these errors will usually be caught downstream. Rather than
utilizing convenience functions in the Python C-API, raw strings are created
and data is copied in using memcpy. This is signficantly faster.


# Comparison

* Biopython [@Cock2009Biopython]
* fastp [@Chen2018Fastp]
* fastq-and-furious [@fastq-and-furious]
* pyfastx [@Du2020Pyfastx]
* pysam [@pysam], wraps htslib [@Bonfield2021HTSlib]
* scikit-bio [@scikit-bio]


# Benchmarks

* Include version numbers for all libraries


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Availability

The homepage for dnaio is at https://dnaio.readthedocs.io/.
It contains installation instructions, a tutorial and the API reference.
The software is  covered by the MIT license.
dnaio is also distributed on Bioconda `[@Gruening2018Bioconda@]`.


# Acknowledgements

MM is financially supported by the Knut and Alice Wallenberg Foundation as part
of the National Bioinformatics Infrastructure Sweden at SciLifeLab.

# References

<!--
A list of key references, including to other software addressing related needs.
Note that the references should include full names of venues, e.g., journals
and conferences, not abbreviations only understood in the context of a specific
discipline.
-->
