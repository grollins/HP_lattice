from numpy import array, zeros, int32, r_, append, sqrt, sum
from .util import vec2coords, check_viability, compute_energy, is_nonsym, \
                  do_shift


DTYPE = int32


class Chain(object):
    """
    The 2D HP lattice chain. The chain is defined by monomers (beads), which can
    be either Hydrophobic (H) or Polar (P). By convention, the first monomer in
    the chain is fixed at the origin on a two-dimensional square lattice.

    Adjacent monomers are connected by bonds. The directionality of the bonds is
    stored in :class:`Vectors`, and each bond-vector is stored as an integer: 
    0 (up), 1 (left), 2 (down) or 3 (right). A chain with :math:`N` monomers has
    :math:`N-1` vectors.
    
    In addition to the vector representation of the chain, :class:`Chain` objects
    also store the 2D coordinates of each monomer. The set of coordinates is stored
    in a :class:`Coords` object.

    :param str hpstring: String that specifies H or P for each monomer.
                         Example: ``PHPPHP``
    :param list initial_vec: List that specifies the direction of each bond.
                             Example: ``[0,0,1,2,1]``, which would correspond
                             to up, up, left, down, left.

    """

    class Vectors(object):
        """
        :class:`Vectors` stores the directionality of the bonds between monomers in a
        :class:`Chain`. Each bond-vector is stored as an integer:
        0 (up), 1 (left), 2 (down) or 3 (right). A chain with :math:`N` monomers
        has :math:`N-1` vectors.

        :param list vec_list: List that specifies the direction of each bond.
                              Example: ``[0,0,1,2,1]``, which would correspond
                              to up, up, left, down, left.

        """
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
            """
            :return: copy of this object
            :rtype: :class:`Vectors`
            """
            v = Chain.Vectors()
            v.vec = self.vec.copy()
            return v

        def set(self, vec):
            """
            Change all of the vectors.

            :param vec: array of vector values
            :type vec: :class:`numpy.ndarray`
            """
            self.vec = vec

        def set_idx(self, idx, value):
            """
            Change one of the vectors.

            :param int idx: the index of the vector to change
            :param int value: the value to change the vector to
            """
            self.vec[idx] = value

        def get(self, idx):
            """
            Get one of the vectors.

            :param int idx: the index of a vector 
            """
            return self.vec[idx]

        def grow(self):
            """
            Append another vector to C-terminus of the chain, pointing up.
            """
            self.vec = r_[self.vec, 0]

        def pop(self):
            """
            Remove the last vector from the C-terminus of the chain.
            """
            self.vec = self.vec[:-1]

        def increment(self, idx):
            """
            Rotate a vector clockwise.

            :param int idx: the index of a vector
            """
            self.vec[idx] += 1

        def as_npy_array(self):
            """
            Convert :class:`Vectors` object to :class:`numpy.ndarray`

            :return: array of vector values
            :rtype: :class:`numpy.ndarray`
            """
            return self.vec


    class Coords(object):
        """
        :class:`Coords` is the set of monomer coordinates for a chain.

        :param int num_monomers: The number of monomers in the chain.

        """
        def __init__(self, num_monomers=0):
            if num_monomers > 0:
                self.coords = zeros([num_monomers,2], DTYPE)
            else:
                self.coords = None

        def __len__(self):
            return len(self.coords)

        def __str__(self):
            return str(self.coords)

        def copy(self):
            """
            :return: a copy of this object
            :rtype: :class:`Coords`
            """
            c = Chain.Coords()
            c.coords = self.coords.copy()
            return c

        def set(self, coords, idx=None):
            """
            Change the coordinates for all monomers or one monomer if
            *idx* is defined.

            :param coords: 2D array of coordinates
            :type coords: :class:`numpy.ndarray`
            :param int idx: (optional), the index of a monomer
            """
            if idx is None:
                self.coords = coords
            else:
                self.coords[idx,0] = coords[0]
                self.coords[idx,1] = coords[1]

        def get(self, idx=None):
            """
            Get the coordinates of a monomer.

            :param int idx: the index of a monomer
            """
            if idx is None:
                return self.coords
            else:
                return self.coords[idx,:]

        def vec2coords(self, vec):
            """
            Convert a :class:`Vectors` object into a set of coordinates,
            and replace the currently stored coordinates with the result.

            :param vec: convert these vectors into a set of coordinates
            :type vec: :class:`Vectors`
            """
            # delegate to faster Cython implementation
            self.coords = vec2coords(vec.as_npy_array(), self.coords)

        def as_npy_array(self):
            """
            Convert :class:`Coords` object to :class:`numpy.ndarray`

            :return: array of coordinates
            :rtype: :class:`numpy.ndarray`
            """
            return self.coords

        def is_viable(self):
            """
            Check for steric overlap and chain crossovers, neither of which
            are allowed.

            :return: ``True`` if no problems are found with the chain.
            :rtype: bool
            """
            return check_viability(self.coords)

        def distance_between_pts(self, idx1, idx2):
            """
            Compute the cartesian distance between two monomers.

            :param int idx1: the index of a monomer 
            :param int idx2: the index of a monomer
            """
            return sqrt(sum((self.coords[idx1,:] - self.coords[idx2,:])**2))

        def grow(self):
            """
            Append another monomer to C-terminus of the chain, one space above
            the current C-terminal monomer.
            """
            c = self.coords[-1,:]
            app_coord = array([[c[0], c[1]+1],], dtype=int32)
            self.coords = append(self.coords, app_coord, axis=0)

        def pop(self):
            """
            Remove the last monomer from the C-terminus of the chain.
            """
            self.coords = self.coords[:-1,:]

        def rotate_up_to_right(self, idx):
            """
            Increments the x-coord and decrements the y-coord of a monomer.

            :param int idx: the index of a monomer
            """
            self.coords[idx,0] += 1
            self.coords[idx,1] -= 1

        def rotate_right_to_down(self, idx):
            """
            Decrements the x-coord and decrements the y-coord of a monomer.

            :param int idx: the index of a monomer
            """
            self.coords[idx,0] -= 1
            self.coords[idx,1] -= 1

        def rotate_down_to_left(self, idx):
            """
            Decrements the x-coord and increments the y-coord of a monomer.

            :param int idx: the index of a monomer
            """
            self.coords[idx,0] -= 1
            self.coords[idx,1] += 1


    def __init__(self, hpstring, initial_vec):
        print '\tInitializing Chain.py object...'
        # The HP sequence as a string
        self.hpstring = hpstring
        # the chain length
        self.n = len(self.hpstring)
        # The HP seq in binary rep (H=1 P=0)
        self.hpbinary = self.hpstr2bin()
        
        H_inds = [idx for idx, bead in enumerate(self.hpstring) if bead == 'H']
        self.H_inds = array(H_inds, int32)

        # an (n-1)-dimensional vector representation of the chain
        self.vec = Chain.Vectors(initial_vec)

        # the 2D coordinates of the chain, as a list of tuples 
        self.coords = Chain.Coords(len(self.vec)+1)
        self.vec2coords()

        # Initialize the vec, coords, and viable of any
        #  **proposed** new conformation.
        # Having these variables is convenient for use with
        # Monte Carlo algorithms, e.g.
        self.nextvec = self.vec.copy()
        self.nextcoords = self.coords.copy()

    def __len__(self):
        return self.n

    def __str__(self):
        return "%s\n%s" % (self.vec, self.coords)

    def get_vec_length(self):
        """
        :return: The number of vectors in the chain, which should be one less
                 than the number of monomers.
        :rtype: int
        """
        return len(self.vec)

    def hpstr2bin(self):
        """
        Convert a string of type ``HPHPHPPPHHP`` to a list of 1s and 0s.
        """
        binseq = []
        for i in range(0, len(self.hpstring)):
            if self.hpstring[i] == 'H':
                binseq.append(1)
            else:
                binseq.append(0)
        return binseq

    def vec2coords(self):
        """
        Convert a :class:`Vectors` object into a set of coordinates,
        and replace the currently stored coordinates with the result.

        :param vec: generate coordinates from current vectors
        :type vec: :class:`Vectors`
        """
        # delegate to the implementation in Coords
        self.coords.vec2coords(self.vec)

    def is_viable(self):
        """
        Check for steric overlap and chain crossovers, neither of which
        are allowed.

        :return: ``True`` if no problems are found with the chain.
        :rtype: bool
        """
        # delegate to the implementation in Coords
        return self.coords.is_viable()

    def contactstate(self):
        """
        Find all contacts between pairs of ``H`` monomers that are separated by at
        least three positions along the chain, i.e. will not include
        ``(i, i+1)`` and ``(i, i+2)`` pairs.

        :return: list of ``(idx1,idx2)`` contacts (tuples)
        :rtype: list
        """
        contactstate = []
        for c in range(0, len(self.coords)-1):
            for d in range((c+3), len(self.coords)):
                if self.hpstring[c] == 'H' and self.hpstring[d] == 'H' and \
                self.coords.distance_between_pts(c, d) == 1:
                    contactstate.append((c,d))
        return contactstate 

    def energy(self, epsilon=0.):
        """
        Compute energy of chain, based on hydrophobic contacts.

        :param float epsilon: the energy of one hydrophobic contact.
        :return: the total energy of the chain.
        :rtype: float
        """
        return compute_energy(epsilon, self.coords.as_npy_array(), self.H_inds)

    def grow(self):
        """
        Add a new monomer to C-terminus of the chain.
        """
        # Add a "0" onto the vec chain...
        self.vec.grow()
        # ... update the coords
        self.coords.grow()

    def shift(self):
        """
        Shifts the chain vector to the 'next' list, according to an
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

        :return: ``1`` if no more shifts are possible and ``0`` otherwise.
        :rtype: int
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
        """
        Many of the conformations are related by rotations and reflections.
        We define a "non-symmetric" conformation to have the first
        direction ``0`` and the first turn be a ``1`` (right turn)
        
        :return: ``1`` if the chain is non-symmetric, ``0`` otherwise.
        :rtype: int
        """
        # delegate to the faster Cython implementation
        return is_nonsym( self.vec.as_npy_array() )

    def is_first_vec_one(self):
        """
        :return: ``True`` if the first vector points to the right
        :rtype: bool
        """
        return (self.vec.get(0) == 1)

    def get_coord_array(self):
        """
        :return: the coordinates of the chain
        :rtype: :class:`numpy.ndarray`
        """
        return self.coords.as_npy_array()

    def get_hp_string(self):
        """
        :return: example ``HPHPHPPPHHP``
        :rtype: string
        """
        return self.hpstring

    def do_three_bead_flip(self, vecindex):
        """
        Swap the values of a pair of neighboring chain vectors.

        :param int vecindex: swap ``vecindex`` with ``vecindex+1``
        """
        # do flip if vec directions are different
        tmp1 = self.nextvec.get(vecindex) 
        tmp2 = self.nextvec.get(vecindex + 1)
        if tmp1 != tmp2:
            self.nextvec.set_idx(vecindex, tmp2)
            self.nextvec.set_idx(vecindex + 1, tmp1)

    def do_crankshaft(self, vecindex):
        """
        Swap the values of a pair of chain vectors that are next-nearest neighbors.

        :param int vecindex: swap ``vecindex`` with ``vecindex+2``
        """
        # do crankshaft if vec directions are different
        tmp1 = self.nextvec.get(vecindex)
        tmp2 = self.nextvec.get(vecindex + 2)
        if tmp1 != tmp2:
            self.nextvec.set_idx(vecindex, tmp2)
            self.nextvec.set_idx(vecindex + 2, tmp1)

    def do_rigid_rot(self, vecindex, direction):
        """
        Rotate a chain vector clockwise or counterclockwise. The possible
        clockwise rotations are:

            0 --> 1

            1 --> 2
        
            2 --> 3
        
            3 --> 0

        :param int vecindex: a valid vector index
        :param int direction: ``1`` for clockwise, ``-1`` for counterclockwise
        """
        self.nextvec.vec[vecindex:] = \
            (self.nextvec.vec[vecindex:] + direction) % 4

    def update_chain(self):
        """
        Accept recent chain move. This is usually called after a trial monte
        carlo move to accept the chain perturbation.
        """
        self.vec.vec[:] = self.nextvec.vec[:]
        self.coords.coords[:,:] = self.nextcoords.coords[:,:]

    def nextviable(self):
        """
        Check for steric overlap and chain crossovers in the not-yet-accepted
        next configuration of the chain. This is usually called before accepting
        a monte carlo move to make sure that the chain doesn't end up in a
        disallowed state.

        :return: ``True`` if no problems are found with the chain.
        :rtype: bool
        """
        self.nextcoords.vec2coords(self.nextvec)
        return self.nextcoords.is_viable()
