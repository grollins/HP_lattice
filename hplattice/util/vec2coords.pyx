import numpy as N
cimport numpy as N
import cython

DTYPE = N.int32
ctypedef N.int32_t DTYPE_t

# Cython compiler directives set for efficiency:
# - No bound checks on index operations
# - No support for negative indices
# - Division uses C semantics
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef N.ndarray[DTYPE_t, ndim=2] vec2coords(N.ndarray[DTYPE_t, ndim=1] chain_vecs, N.ndarray[DTYPE_t, ndim=2] coords):
    """Convert an array of chain vectors to a 2d array of
       coordinates."""
    cdef int cv_max = chain_vecs.shape[0]
    cdef int x, y, i
    x = 0
    y = 0
    coords[0,0] = x
    coords[0,1] = y
    for i in range(cv_max):
        if chain_vecs[i] == 0:
            y = y + 1
        if chain_vecs[i] == 1:
            x = x + 1
        if chain_vecs[i] == 2:
            y = y - 1
        if chain_vecs[i] == 3:
            x = x - 1
        coords[i+1,0] = x
        coords[i+1,1] = y
    return coords
