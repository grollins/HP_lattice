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
cpdef N.ndarray[DTYPE_t, ndim=2] shift(N.ndarray[DTYPE_t, ndim=1] vec, N.ndarray[DTYPE_t, ndim=2] coords):
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

    vec = self.vec.tolist()
    coords = [(row[0], row[1]) for row in self.coords]

    # pop off all the trailing 3's      
    while 1:
        if len(vec) > 0 and vec[-1] == 3:
            vec.pop() # update vec
            coords.pop() # update coords
        else:
            break

    ### update viability
    self.viable = 1
    for c in coords:
        if coords.count(c) > 1:
            self.viable = 0

    ## if popping of all the '3's left an empty vec list, then we're done!
    if len(vec) == 0:
        self.viable = 0
        self.vec = array(vec, int32)
        self.coords = self.vec2coords(self.vec, self.coords)
        return(1)

    i = len(vec) - 1 # the last vec index

    # Otherwise, increment the remaining non-"3" elements
    if vec[i] == 0:
        # update vec 
        vec[i] = vec[i] + 1
        # update coords: rotate "up" to "right"
        coords[i+1] = (coords[i+1][0] + 1, coords[i+1][1] - 1)
        # update viability
        self.viable = 1
        if coords.count(coords[i+1]) > 1:
            self.viable = 0

    elif vec[i] == 1:
        # update vec
        vec[i] = vec[i] + 1
        # update coords: rotate "right" to "down"
        coords[i+1] = (coords[i+1][0] - 1, coords[i+1][1] - 1)
        # update viability
        self.viable = 1
        if coords.count(coords[i+1]) > 1:
            self.viable = 0

    else:
        # i.e. self.vec[i] == 2:
        # update vec
        vec[i] = vec[i] + 1
        # update coords: rotate "down" to "left"
        coords[i+1] = (coords[i+1][0] - 1, coords[i+1][1] + 1)
        # update viability
        self.viable = 1
        if coords.count(coords[i+1]) > 1:
            self.viable = 0
    
    self.vec = array(vec, int32)
    self.coords = self.vec2coords(self.vec, self.coords)
    return(0)
