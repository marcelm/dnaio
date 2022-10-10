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

- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.


# Summary

<!--
A summary describing the high-level functionality and purpose of the
software for a diverse, non-specialist audience.
-->

Low-cost, high-throughput DNA sequencing has become the workhorse of life
science research. Compressed FASTQ files [@Cock2010Sanger] are the nearly
universally used file format for representing the nucleotide sequences
obtained from a sequencing machine.

...

# Statement of need

<!--
A Statement of need section that clearly illustrates the research
purpose of the software and places it in the context of related work.
-->


# To Do

* already used by Cutadapt and therefore has >1000 daily downloads
* mention license
* homepage with tutorial and API reference at https://dnaio.readthedocs.io/
* test coverage
* benchmarks and comparison to other tools
* features compared to other tools (feature matrix)
* tricks used to make FASTQ parsing fast


# Alternatives

* Biopython [@Cock2009Biopython]
* fastp [@Chen2018Fastp]
* fastq-and-furious [@fastq-and-furious]
* pyfastx [@Du2020Pyfastx]
* pysam [@pysam], wraps htslib [@Bonfield2021HTSlib]
* scikit-bio [@scikit-bio]



# Benchmarks

* Include version numbers for all libraries


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

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
