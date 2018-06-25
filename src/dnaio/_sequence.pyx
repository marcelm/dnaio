# cython: profile=False, language_level=3, emit_code_comments=False
from ._util import shorten
from .exceptions import FileFormatError


cdef class Sequence(object):
	"""
	A record in a FASTA or FASTQ file. For FASTA, the qualities attribute
	is None. For FASTQ, qualities is a string and it contains the qualities
	encoded as ascii(qual+33).
	"""

	def __cinit__(self, str name, str sequence, str qualities=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities

		if qualities is not None and len(qualities) != len(sequence):
			rname = shorten(name)
			raise FileFormatError("In read named {0!r}: length of quality sequence ({1}) and length "
				"of read ({2}) do not match".format(
					rname, len(qualities), len(sequence)))

	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(shorten(self.qualities))
		return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(shorten(self.name), shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __richcmp__(self, other, int op):
		if 2 <= op <= 3:
			eq = self.name == other.name and \
				self.sequence == other.sequence and \
				self.qualities == other.qualities
			if op == 2:
				return eq
			else:
				return not eq
		else:
			raise NotImplementedError()

	def __reduce__(self):
		return (Sequence, (self.name, self.sequence, self.qualities))
