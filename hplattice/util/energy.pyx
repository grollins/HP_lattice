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
cpdef double energy(double epsilon, N.ndarray[DTYPE_t, ndim=2] coords, N.ndarray[DTYPE_t, ndim=1] H_inds):
    """Calculate potential energy of the chain."""
    cdef int i, j, d
    cdef double E
    cdef int nc = 0
    cdef int h_max = H_inds.shape[0]

    for i in range(h_max):
        for j in range(i+1, h_max):
            if (H_inds[j] - H_inds[i]) >= 3:
                d = (coords[i,0] - coords[j,0])**2 + (coords[i,1] - coords[j,1])**2
                if d == 1:
                    nc += 1
    E = nc * epsilon
    return E
