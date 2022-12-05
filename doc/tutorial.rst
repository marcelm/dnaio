Tutorial
========

This should get you started with using ``dnaio``.
The only essential concepts to know about are
the `dnaio.open` function and the `~dnaio.SequenceRecord` object.


Reading
-------

The main interface for reading and writing sequence files is the `dnaio.open` function.
For example, this program reads in a FASTQ file and computes the total number of nucleotides
it contains::

    import dnaio

    with dnaio.open("reads.fastq.gz") as reader:
        bp = 0
        for record in reader:
            bp += len(record)
    print(f"The input file contains {bp/1E6:.1f} Mbp")

As can be seen from the ``.gz`` file extension,
the input file is gzip-compressed.
`dnaio.open` detects and handles this automatically by opening the file with
`xopen <https://github.com/pycompression/xopen/>`_.

Here, the call to `dnaio.open` returns a `~dnaio.FastqReader` object.
Iterating over it in the ``for`` loop results in `~dnaio.SequenceRecord` objects.
Calling ``len()`` on a ``SequenceRecord`` returns the number of
nucleotides in the record.

A ``SequenceRecord`` has the attributes ``name``, ``sequence``
and ``qualities``. All of these are ``str`` objects.
The ``qualities`` attribute is ``None`` when reading FASTA files.

The following program uses the ``name`` attribute
to check whether any sequence names are duplicated in a FASTA file::

    import dnaio

    seen = set()
    with dnaio.open("sequences.fasta") as reader:
        for record in reader:
            if record.name in seen:
                print(record.name, "is duplicated")
            seen.add(record.name)

Writing
-------

To open a sequence file for writing,
pass the ``mode="w"`` argument to ``dnaio.open``::

    import dnaio

    with dnaio.open("onerecord.fastq.gz", mode="w") as writer:
        writer.write(dnaio.SequenceRecord("name", "ACGT", "#B!#"))

Here, a `~dnaio.FastqWriter` object is returned by ``dnaio.open``,
which has a `~dnaio.FastqWriter.write()` method that accepts a ``SequenceRecord``.

A possibly more common use case is to read an input file,
modify the reads and write them to a new output file.
The following example program shows how that can be done.
It truncates all reads in the input file to a length of 30 nt
and writes them to another file::

    import dnaio

    with dnaio.open("in.fastq.gz") as reader, dnaio.open("out.fastq.gz", mode="w") as writer:
        for record in reader:
            record = record[:30]
            writer.write(record)

This also shows that `~dnaio.SequenceRecord` objects support slicing:
``record[:30]`` returns a new ``SequenceRecord`` object with the sequence and qualities
trimmed to the first 30 characters, leaving the name unchanged.


Paired-end data
---------------

Paired-end data is supported in two forms: Two separate files or interleaved.

To read from separate files, provide two input file names to the ``dnaio.open``
function::

    import dnaio

    with dnaio.open("reads.1.fastq.gz", "reads.2.fastq.gz") as reader:
        bp = 0
        for r1, r2 in reader:
            bp += len(r1) + len(r2)
        print(f"The paired-end input contains {bp/1E6:.1f} Mbp")

Here, ``dnaio.open`` returns a `~dnaio.TwoFilePairedEndReader`.
It also supports iteration, but instead of a plain ``SequenceRecord``,
it returns a tuple of two ``SequenceRecord`` instances.

To read from interleaved paired-end data,
pass ``interleaved=True`` to ``dnaio.open`` instead of a second file name::

    ...
    with dnaio.open("reads.interleaved.fastq.gz", interleaved=True) as reader:
        ...

The ``PairedEndReader`` classes check whether the input files are properly paired,
that is, whether they have the same number of reads in both inputs and whether the
read names match.
For this reason, always use a single call to ``dnaio.open`` to open paired-end files
(that is, avoid opening them as two single-end files.)

To demonstrate how to write paired-end data,
we show a program that reads from a single-end FASTQ file and converts the records to
simulated paired-end reads by writing the first 30 nt to R1 and the last 30 nt
to R2::

    import dnaio

    with dnaio.open("in.fastq.gz") as reader, \
            dnaio.open("out.1.fastq.gz", "out.2.fastq.gz", mode="w") as writer:
        for record in reader:
            r1 = record[:30]
            r2 = record[-30:]
            writer.write(r1, r2)

The ``writer`` in this case is a `~dnaio.TwoFilePairedEndWriter`
and its `~dnaio.TwoFilePairedEndWriter.write()` method
expects two ``SequenceRecord`` arguments.
