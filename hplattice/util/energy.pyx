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
cpdef energy(double epsilon, N.ndarray[DTYPE_t, ndim=2] coords, N.ndarray[DTYPE_t, ndim=1] H_inds):
    """Calculate potential energy of the chain."""
    cdef int i, j, d, h1, h2
    cdef double E
    cdef int nc = 0
    cdef int h_max = H_inds.shape[0]
    contacts = []

    for i in range(h_max):
        for j in range(i+1, h_max):
            h2 = H_inds[j]
            h1 = H_inds[i]
            if (h2 - h1) >= 3:
                d = (coords[h1,0] - coords[h2,0])**2 + (coords[h1,1] - coords[h2,1])**2
                if d == 1:
                    contacts.append((h1,h2))
                    nc += 1
    E = nc * epsilon
    return E, contacts
