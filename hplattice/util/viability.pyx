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
cpdef int viability(N.ndarray[DTYPE_t, ndim=2] coords):
    """Return 1 if the chain coordinates are self-avoiding,
           0 if not."""
    cdef int i, j, d
    cdef int row_max = coords.shape[0]
    for i in range(row_max):
        for j in range(i+1, row_max):
            d = (coords[i,0] - coords[j,0])**2 + (coords[i,1] - coords[j,1])**2
            if d < 1.0:
                return 0
    return 1
