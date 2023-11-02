The dnaio API
=============


.. module:: dnaio


The open function
-----------------

.. autofunction:: open


The ``SequenceRecord`` class
------------------------------

.. autoclass:: dnaio.SequenceRecord
   :members:
   :special-members: __len__, __getitem__

   .. automethod:: __init__(name: str, sequence: str, qualities: Optional[str] = None)


Reader and writer interfaces
----------------------------

.. autoclass:: SingleEndReader
   :members: __iter__

.. autoclass:: PairedEndReader
   :members: __iter__

.. autoclass:: SingleEndWriter
   :members: write

.. autoclass:: PairedEndWriter
   :members: write

.. autoclass:: MultipleFileWriter
   :members: write, write_iterable


Reader and writer classes
-------------------------

The `dnaio.open` function returns an instance of one of the following classes.
They can also be used directly if needed.


.. autoclass:: FastaReader
   :show-inheritance:

.. autoclass:: FastaWriter
   :show-inheritance:

.. autoclass:: FastqReader
   :show-inheritance:

.. autoclass:: FastqWriter
   :show-inheritance:

.. autoclass:: BamReader
   :show-inheritance:

.. autoclass:: TwoFilePairedEndReader
   :show-inheritance:

.. autoclass:: TwoFilePairedEndWriter
   :show-inheritance:

.. autoclass:: InterleavedPairedEndReader
   :show-inheritance:

.. autoclass:: InterleavedPairedEndWriter
   :show-inheritance:

.. autoclass:: MultipleFileReader
   :members: __iter__

.. autoclass:: MultipleFastaWriter
   :show-inheritance:

.. autoclass:: MultipleFastqWriter
   :show-inheritance:


Chunked reading of sequence records
-----------------------------------

The following functions can be used to very quickly split up the input file(s)
into similarly-sized chunks without actually parsing the records. The chunks
can then be distributed to worker threads or subprocesses and be parsed and
processed there.

.. autofunction:: read_chunks
.. autofunction:: read_paired_chunks


Functions
---------

.. autofunction:: records_are_mates


Exceptions
----------

.. autoexception:: UnknownFileFormat

.. autoexception:: FileFormatError

.. autoexception:: FastaFormatError
   :show-inheritance:

.. autoexception:: FastqFormatError
   :show-inheritance:
