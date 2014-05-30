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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int is_nonsym(N.ndarray[DTYPE_t, ndim=1] vec):
    """Many of the conformations are related by rotations and reflections.
       We define a "non-symmetric" conformation to have the first
       direction '0' and the first turn be a '1' (right turn)
       nonsym() returns 1 if the vec list is non-symmetric, 0 otherwise
    """
    cdef int i, v
    cdef int vmax = vec.shape[0]
    if vmax > 0:
        # walk along chain until you get to the first non '0' vec:
        # for i in range(vmax - 1):
        #     if vec[i] > 0:
        #         break
        i = 0
        v = vec[0]
        while (v == 0) and (i < (vmax - 1)):
            i = i + 1
            v = vec[i]
        if vec[0] == 0:
            if (vec[i] == 1) or (vec[i] == 0):
                return 1
    return 0

    # if len(self.vec) > 0:
    #     # walk along chain until you get to the first non '0' vec:
    #     i = 0
    #     v = self.vec.get(0)
    #     while (v == 0) and (i < (len(self.vec) - 1)):
    #         i = i + 1
    #         v = self.vec.get(i)
    #     if self.vec.get(0) == 0:
    #         if (v == 1) or (v == 0):
    #             return 1
    # return 0