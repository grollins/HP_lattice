from random import random
from math import floor, exp
from numpy import array, int32


BOLTZ_CONST = 0.001987  # (kcal/K.mol) Boltzmann's constant


class Monty(object):
    """A collection of functions to perform Monte Carlo move-set operations
       on an HP lattice Chain object."""

    def __init__(self, config, temp, chain):
        """Initialize the Monte Carlo object..."""
        print '\tcreating Monty.py object....'
        # indices of hydrophic beads
        H_inds = [idx for idx, bead in enumerate(chain.hpstring) if bead == 'H']
        self.H_inds = array(H_inds, int32)

        # The names of the available Monte Carlo movesets
        self.movesets = ['MC1','MC2','MC3','MC4']
        # The temperature (in K)
        self.temp = temp
        # The replica number with this temperature
        self.tempfromrep = config.REPLICATEMPS.index(temp)
        # The energetic strength of a contact
        self.epsilon = config.epsilon
        # Boltzmann's constant
        # self.k = config.k

        # (both of these are copies from Config() )
        self.restraint = DistRestraint(config.RESTRAINED_STATE, config.KSPRING)
        self.lastenergy = self.energy(chain) + self.restraint.energy(chain)

    def kT(self):
        return BOLTZ_CONST * self.temp

    def move1(self, chain, vecindex=None, direction=None):
        """Apply moveset 'MC1' to the chain:
           (i)  three-bead flips
           (ii) end flips
           REFERENCE: Dill and Chan, 1994, 1996.
        """
        if vecindex is None:
            r = random()
            vecindex = int(floor((chain.n - 1.0001)*r))
        else:
            pass

        if direction is None:
            s = random()
            if s < 0.5:
                direction = 1
            else:
                direction = -1
        else:
            pass
 
        if (vecindex == 0) or (vecindex == len(chain.nextvec)-1):
            # if the vec index is on the end, do an end flip
            chain.do_rigid_rot(vecindex, direction)

        else:
            # do a three-bead flip, i.e. switch two adjacent {n,e,w,s} directions
            chain.do_three_bead_flip(vecindex)

    def move2(self, chain, vecindex=None, direction=None, moveseed=None):
        """
        Apply moveset MC2 to the chain:
        (i)   three-bead flips
        (ii)  end flips
        (iii) crankshaft moves
        (iv)  rigid rotations
        REFERENCE:  Dill and Chan, 1994, 1996
        """
        if vecindex is None:
            r = random()
            vecindex = int(floor((chain.n - 1.0001)*r))
        else:
            pass

        if direction is None:
            s = random()
            if s < 0.5:
                direction = 1
            else:
                direction = -1
        else:
            pass

        length_of_nextvec = chain.get_vec_length()
        
        if moveseed is None:
            t = random()
        else:
            t = moveseed

        # 1/3 of the time do a three-bead flip
        if (t < 0.33333) and (vecindex < length_of_nextvec-1):
            flip_successful = chain.do_three_bead_flip(vecindex)
            if flip_successful:
                pass
            else:
                chain.do_rigid_rot(vecindex, direction)

        # 1/3 of the time do a crankshaft
        elif (t < 0.66666) and (vecindex < length_of_nextvec-2):
            crank_successful = chain.do_crankshaft(vecindex)
            if crank_successful:
                pass
            else: 
                chain.do_rigid_rot(vecindex, direction)

        # 1/3 of the time do a rigid rotation
        else: 
            chain.do_rigid_rot(vecindex, direction)

    def move3(self, chain, vecindex=None, direction=None):
        """
        Apply moveset 'MC3' to the chain.
        This is just a simple set to change the direction of
        a single chain link.
        Example:
            [0,0,0,0,0] --> [0,0,1,0,0]
        where {0,1,2,3}={n,e,s,w} direction
        About 5% viable moves are expected.
        """
    
        if vecindex is None:
            r = random()
            vecindex = int(floor((chain.n - 1.0001)*r))
        else:
            pass

        if direction is None:
            s = random()
            if s < 0.5:
                direction = 1
            else:
                direction = -1
        else:
            pass
  
        chain.do_rigid_rot(vecindex, direction)

    def metropolis(self, replica):
        """
        Accept Chain.nextvec over Chain.vec according to a Metropolis
        criterion.
        """
        randnum = random()

        # accept with Metroplis criterion
        thisenergy = self.energy(replica.chain) + \
                        self.restraint.energy(replica.chain)  
        boltzfactor = exp( (thisenergy - self.lastenergy) / self.kT() )

        if randnum < boltzfactor:
            # update the chain
            replica.chain.update_chain()
            # update the lastenergy
            self.lastenergy = thisenergy
            return 1
        else:
            return 0

    def energy(self, chain):
        E, state = chain.energy(self.epsilon)
        return E


class DistRestraint:
    """ For now, this is a harmonic constraint over a squared distance D = d^2
     where D = sum_{i,j} d^2_ij over all contacts."""
    def __init__(self, contacts, kspring):
        """Initialize the DistRestraint object"""
        # a list of tuples (start of chain is 0)
        self.contacts = contacts
        # spring constant for restraint
        # (J/[lattice space]^2)
        self.kspring = kspring

    def energy(self, chain):
        """ return the energy of the distance restraint"""
        return self.kspring * self.D(chain)

    def D(self, chain):
        """Return the sum of squared-distances over the selected contacts."""
        D = 0.0
        coords = chain.get_coord_array()
        for i in range(0, len(self.contacts)):
            # print 'contact', i, self.contacts[i], chain.coords
            c = self.contacts[i][0]
            d = self.contacts[i][1]
            D = D + (coords[c][0]-coords[d][0])**2
            D = D + (coords[c][1]-coords[d][1])**2
        return D
