from numpy import array, zeros, int32, nonzero, r_, append, sqrt, sum
from .util import vec2coords, check_viability, compute_energy, is_nonsym, \
                  do_shift


DTYPE = int32


class Chain:
    """
    An object to represent the 2D HP lattice chain and its attributes,
    with method functions.
    """

    class Vectors(object):
        """docstring for Vectors"""
        def __init__(self, vec_list=None):
            if vec_list:
                self.vec = array(vec_list, DTYPE)
            else:
                self.vec = None

        def __len__(self):
            return len(self.vec)

        def __str__(self):
            return str(self.vec)

        def copy(self):
            v = Chain.Vectors()
            v.vec = self.vec.copy()
            return v

        def set(self, vec):
            self.vec = vec

        def set_idx(self, idx, value):
            self.vec[idx] = value

        def get(self, idx):
            return self.vec[idx]

        def grow(self):
            self.vec = r_[self.vec, 0]

        def pop(self):
            self.vec = self.vec[:-1]

        def increment(self, idx):
            self.vec[idx] += 1

        def as_npy_array(self):
            return self.vec


    class Coords(object):
        """docstring for Coords"""
        def __init__(self, num_points=0):
            if num_points > 0:
                self.coords = zeros([num_points,2], DTYPE)
            else:
                self.coords = None

        def __len__(self):
            return len(self.coords)

        def __str__(self):
            return str(self.coords)

        def copy(self):
            c = Chain.Coords()
            c.coords = self.coords.copy()
            return c

        def set(self, coords):
            self.coords = coords

        def vec2coords(self, vec):
            self.coords = vec2coords(vec.as_npy_array(), self.coords)

        def as_npy_array(self):
            return self.coords

        def is_viable(self):
            return check_viability(self.coords)

        def distance_between_pts(self, idx1, idx2):
            return sqrt(sum((self.coords[idx1,:] - self.coords[idx2,:])**2))

        def grow(self):
            c = self.coords[-1,:]
            app_coord = array([[c[0], c[1]+1],], dtype=int32)
            self.coords = append(self.coords, app_coord, axis=0)

        def pop(self):
            self.coords = self.coords[:-1,:]

        def rotate_up_to_right(self, idx):
            self.coords[idx,0] += 1
            self.coords[idx,1] -= 1

        def rotate_right_to_down(self, idx):
            self.coords[idx,0] -= 1
            self.coords[idx,1] -= 1

        def rotate_down_to_left(self, idx):
            self.coords[idx,0] -= 1
            self.coords[idx,1] += 1


    def __init__(self, config):
        print '\tInitializing Chain.py object...'
        # The HP sequence as a string
        self.hpstring = config.HPSTRING.strip()
        # the chain length
        self.n = len(self.hpstring)
        # The HP seq in binary rep (H=1 P=0)
        self.hpbinary = self.hpstr2bin()
        
        H_inds = [idx for idx, bead in enumerate(self.hpstring) if bead == 'H']
        self.H_inds = array(H_inds, int32)

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
        # self.vec = array(config.INITIALVEC, int32)
        self.vec = Chain.Vectors(config.INITIALVEC)

        # the 2D coordinates of the chain, as a list of duples 
        # self.coords = zeros([len(self.vec)+1,2], int32)
        self.coords = Chain.Coords(len(self.vec)+1)
        self.vec2coords()

        # Initialize the vec, coords, and viable of any
        #  **proposed** new conformation.
        # Having these variables is convenient for use with
        # Monte Carlo algorithms, e.g.
        self.nextvec = self.vec.copy()
        self.nextcoords = self.coords.copy()

    def __str__(self):
        return "%s\n%s" % (self.vec, self.coords)

    def get_vec_length(self):
        return len(self.vec)

    def hpstr2bin(self):
        """Convert a string of type 'HPHPHPPPHHP' to a list of 1s and 0s.""" 
        binseq = []
        for i in range(0, len(self.hpstring)):
            if self.hpstring[i] == 'H':
                binseq.append(1)
            else:
                binseq.append(0)
        return binseq

    def vec2coords(self):
        self.coords.vec2coords(self.vec)

    def is_viable(self):
        return self.coords.is_viable()

    def contactstate(self):
        """Return the contact state of the chain as a list of
           (res1,res2) contacts (tuples),
           where the residue numbering starts at 0."""
        contactstate = []
        for c in range(0, len(self.coords)-1):
            for d in range((c+3), len(self.coords)):
                if self.hpstring[c] == 'H' and self.hpstring[d] == 'H' and \
                self.coords.distance_between_pts(c, d) == 1:
                    contactstate.append((c,d))
        return contactstate 

    def energy(self, epsilon=0.):
        return compute_energy(epsilon, self.coords.as_npy_array(), self.H_inds)

    def grow(self):
        """
        Add a new link onto the chain vector,
        updating the coords and viability correspondingly.
        """
        # Add a "0" onto the vec chain...
        self.vec.grow()
        # ... update the coords
        self.coords.grow()

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
        
        # pop off all the trailing 3's		
        while 1:
            if len(self.vec) > 0 and self.vec.get(-1) == 3:
                self.vec.pop() # update vec
                self.coords.pop() # update coords
            else:
                break

        is_done, vec, coords = \
            do_shift(self.vec.as_npy_array(), self.coords.as_npy_array())
        self.vec.set(vec)
        self.coords.set(coords)
        return is_done

    def nonsym(self):
        """Many of the conformations are related by rotations and reflections.
           We define a "non-symmetric" conformation to have the first
           direction '0' and the first turn be a '1' (right turn)
           nonsym() returns 1 if the vec list is non-symmetric, 0 otherwise
        """
        return is_nonsym( self.vec.as_npy_array() )

    def is_first_vec_one(self):
        return (self.vec.get(0) == 1)

    def get_coord_array(self):
        return self.coords.as_npy_array()

    def get_hp_string(self):
        return self.hpstring

    def do_three_bead_flip(self, vecindex):
        # do flip if vec directions are different
        tmp1 = self.nextvec.get(vecindex) 
        tmp2 = self.nextvec.get(vecindex+1)
        if tmp1 != tmp2:
            self.nextvec.set_idx(vecindex, tmp2)
            self.nextvec.set_idx(vecindex + 1, tmp1)

    def do_crankshaft(self, vecindex):
        # do crankshaft if vec directions are different
        tmp1 = self.nextvec.get(vecindex)
        tmp2 = self.nextvec.get(vecindex+2)
        if tmp1 != tmp2:
            self.nextvec.set_idx(vecindex, tmp2)
            self.nextvec.set_idx(vecindex + 2, tmp1)

    def do_rigid_rot(self, vecindex, direction):
        self.nextvec.vec[vecindex:] = \
            (self.nextvec.vec[vecindex:] + direction) % 4

    def update_chain(self):
        self.vec.vec[:] = self.nextvec.vec[:]
        self.coords.coords[:,:] = self.nextcoords.coords[:,:]

    def nextviable(self):
        self.nextcoords.vec2coords(self.nextvec)
        return self.nextcoords.is_viable()
