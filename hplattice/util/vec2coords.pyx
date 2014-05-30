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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef shift(N.ndarray[DTYPE_t, ndim=1] vec, N.ndarray[DTYPE_t, ndim=2] coords):
    """Shifts the chain vector to the 'next' list, according to an
       enumeration scheme where the most distal chain vector is
       incremented 0->1->2->3.  After 3, the most distal vector
       element is removed, and the next most distal element is
       incremented.  If there are multiple "3" vectors, this process
       is done recursively.

       Example:
        [0,0,0,0] --> [0,0,0,1] 
        [0,0,1,2] --> [0,0,1,3] 
        [0,1,0,3] --> [0,1,1] 
        [0,3,3,3] --> [1]

        This operation is very useful for enumerating the full space
        of chain conformations.
        shift()  will also update the coords and the viability, accordingly.

       RETURN VALUES
        returns 1 if its the last possible "shift" --> i.e. if it's all
        3's, the search is done
        returns 0 otherwise
    """
    cdef int i
    cdef int vmax = vec.shape[0]

    # while 1:
    #     print 'a', vec, vec[vmax-1]
    #     if vmax > 0 and vec[vmax-1] == 3:
    #         vec = vec[:vmax-1]
    #         coords = coords[:vmax-1,:]
    #         vmax -= 1
    #         print 'b', vec, vec[vmax-1]
    #     else:
    #         break

    ## if popping of all the '3's left an empty vec list, then we're done!
    if vmax == 0:
        return 1, vec, coords

    i = vmax - 1 # the last vec index

    # Otherwise, increment the remaining non-"3" elements
    if vec[i] == 0:
        # update vec 
        vec[i] += 1
        # update coords: rotate "up" to "right"
        coords[i+1,0] += 1
        coords[i+1,1] -= 1

    elif vec[i] == 1:
        # update vec
        vec[i] += 1
        # update coords: rotate "right" to "down"
        coords[i+1,0] -= 1
        coords[i+1,1] -= 1

    else:
        # i.e. self.vec[i] == 2:
        # update vec
        vec[i] += 1
        # update coords: rotate "down" to "left"
        coords[i+1,0] -= 1
        coords[i+1,1] += 1

    return 0, vec, coords
