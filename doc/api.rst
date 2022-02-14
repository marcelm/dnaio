dnaio API
=========

.. module:: dnaio

The open function
-----------------

.. autofunction:: dnaio.open


The ``Sequence`` classes
------------------------

.. autoclass:: dnaio.Sequence
   :members:
   :special-members: __len__, __getitem__

   .. automethod:: __init__(name: str, sequence: str, qualities: Optional[str] = None)

.. autoclass:: dnaio.BytesSequence
   :members:
   :special-members: __len__, __getitem__

   .. automethod:: __init__(name: bytes, sequence: bytes, qualities: Optional[bytes] = None)

Exceptions
----------

.. autoexception:: UnknownFileFormat

.. autoexception:: FileFormatError

.. autoexception:: FastaFormatError

.. autoexception:: FastqFormatError


Reader and writer classes
-------------------------

.. autoclass:: FastaReader

.. autoclass:: FastaWriter
   :members: write

.. autoclass:: FastqReader

.. autoclass:: FastqWriter
   :members: writeseq

   .. py:method:: write(record: Sequence) -> None:

      Write a Sequence record to the FASTQ file.

.. autoclass:: InterleavedPairedEndReader

.. autoclass:: InterleavedPairedEndWriter
   :members: write

.. autoclass:: SingleEndReader

.. autoclass:: PairedEndReader

.. autoclass:: SingleEndWriter

.. autoclass:: PairedEndWriter
   :members: write


Functions
---------

.. autofunction:: read_chunks
.. autofunction:: read_paired_chunks
.. autofunction:: record_names_match
