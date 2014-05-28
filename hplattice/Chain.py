from numpy import array, zeros, int32, nonzero, r_, append
from .util import vec2coords, check_viability


class Chain:
    """
    An object to represent the 2D HP lattice chain and its attributes,
    with method functions.
    """ 

    def __init__(self, config):
        print '\tInitializing Chain.py object...'
        # The HP sequence as a string
        self.hpstring = config.HPSTRING.strip()
        # the chain length
        self.n = len(self.hpstring)
        # The HP seq in binary rep (H=1 P=0)
        self.hpbinary = self.hpstr2bin()
        
        # The chain is represented by an (n-1)-dimensional vector,
        #  which is stored in a python list.
        # Each element represents the direction of successive links
        #  (bonds) along the chain length:
        #
        #     0     up
        #     1     right
        #     2     down
        #     3     left
        #
        # By convention, the first monomer (bead) in the chain is fixed at
        #  the origin on a two-dimensional square lattice.

        # an (n-1)-dimensional vector representation of the chain
        self.vec = array(config.INITIALVEC, int32)

        # the 2D coordinates of the chain, as a list of duples 
        self.coords = zeros([len(self.vec)+1,2], int32)
        self.coords = self.vec2coords(self.vec, self.coords)
        # viable chains are ones which are self-avoiding
        # 1 if viable, 0 if not       
        self.viable = self.viability(self.coords)     	
        
        # Initialize the vec, coords, and viable of any
        #  **proposed** new conformation.
        # Having these variables is convenient for use with
        # Monte Carlo algorithms, e.g.
        self.nextvec = self.vec.copy()
        self.nextcoords = self.coords.copy()
        self.nextviable = self.viable

    def hpstr2bin(self):
        """Convert a string of type 'HPHPHPPPHHP' to a list of 1s and 0s.""" 
        binseq = []
        for i in range(0, len(self.hpstring)):
            if self.hpstring[i] == 'H':
                binseq.append(1)
            else:
                binseq.append(0)
        return binseq

    def vec2coords(self, vec, coords):
        c = vec2coords(vec, coords)
        return c

    def viability(self, thesecoords):
        # """Return 1 if the chain coordinates are self-avoiding,
        #    0 if not."""
        # self.viable = 1
        # coords = [(row[0], row[1]) for row in thesecoords]
        # for c in coords:
        #     if coords.count(c) > 1:
        #         self.viable = 0
        #         break
        self.viable = check_viability(thesecoords)
        return self.viable

    def contactstate(self):
        """Return the contact state of the chain as a list of
           (res1,res2) contacts (tuples),
           where the residue numbering starts at 0."""
        self.coords = self.vec2coords(self.vec, self.coords)
        contactstate = []
        for c in range(0, len(self.coords)-1):
            for d in range((c+3), len(self.coords)):
                if self.hpstring[c] == 'H':
                    if self.hpstring[d] == 'H':
                        if (abs(self.coords[c][0]-self.coords[d][0]) + \
                            abs(self.coords[c][1]-self.coords[d][1])) == 1:
                            contactstate.append((c,d))
        return contactstate 

    def grow(self):
        """
        Add a new link onto the chain vector,
        updating the coords and viability correspondingly.
        """
        # Add a "0" onto the vec chain...
        self.vec = r_[self.vec, 0]
        # ... update the coords
        i = len(self.coords) - 1
        last_coord = self.coords[i,:]
        app_coord = array([[last_coord[0], last_coord[1]+1],], dtype=int32)
        self.coords = append(self.coords, app_coord, axis=0)
        # ... update the viability
        self.viability(self.coords)

    def shift(self):
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
                self.coords = self.coords[:-1,:]
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

    def nonsym(self):
        """Many of the conformations are related by rotations and reflections.
           We define a "non-symmetric" conformation to have the first
           direction '0' and the first turn be a '1' (right turn)
           nonsym() returns 1 if the vec list is non-symmetric, 0 otherwise
        """
        if len(self.vec) > 0:
            # walk along chain until you get to the first non '0' vec:
            i = 0
            while (self.vec[i] == 0) & (i < len(self.vec) - 1) :
                i = i + 1
            if self.vec[0] == 0:
                if (self.vec[i] == 1) | (self.vec[i] == 0):
                    return(1)
        return(0)
