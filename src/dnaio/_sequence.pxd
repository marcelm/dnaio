# cython: profile=False, language_level=3, emit_code_comments=False

cdef class Sequence:
    cdef:
        public str name
        public str sequence
        public str qualities


