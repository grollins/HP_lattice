from random import random
from math import floor, exp
from numpy import array, int32

from .util import compute_energy


class Monty:
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
        self.k = config.k

        # (both of these are copies from Config() )
        self.restraint = DistRestraint(config.RESTRAINED_STATE, config.KSPRING)
        self.lastenergy = self.energy(chain) + self.restraint.energy(chain)

    def _do_rigid_rot(self, vecindex, direction, chain):
        chain.nextvec[vecindex:] = \
            (chain.nextvec[vecindex:] + direction) % 4

    def move1(self, chain):
        """Apply moveset 'MC1' to the chain:
           (i)  three-bead flips
           (ii) end flips
           REFERENCE: Dill and Chan, 1994, 1996.
        """
        r = random()
        s = random()
        vecindex = int(floor((chain.n - 1.0001)*r))
        if s < 0.5:
            direction = 1
        else:
            direction = -1      
 
        ### if the vec index is on the end, do an end flip
        if (vecindex == 0) | (vecindex == len(chain.nextvec)-1):
            chain.nextvec[vecindex] = (chain.nextvec[vecindex] + direction) % 4
        ### else, do a three-bead flip, i.e. switch two adjacent {n,e,w,s} directions
        else:
            tmp = chain.nextvec[vecindex] 
            chain.nextvec[vecindex] = chain.nextvec[vecindex+1]
            chain.nextvec[vecindex + 1] = tmp

        chain.nextcoords = \
            chain.vec2coords(chain.nextvec, chain.nextcoords)
        chain.nextviable = \
            chain.viability(chain.nextcoords)

    def move2(self, chain):
        """
        Apply moveset MC2 to the chain:
        (i)   three-bead flips
        (ii)  end flips
        (iii) crankshaft moves
        (iv)  rigid rotations
        REFERENCE:  Dill and Chan, 1994, 1996
        """
        r = random()
        s = random()
        t = random()    
        vecindex = int(floor((chain.n - 1.0001)*r))
        
        if s < 0.5:
            direction = 1      
        else:
            direction = -1      

        # if possible, 1/3 of the time do a three-bead flip
        # (dirs must be different)
        if (t < 0.33333) and (vecindex < len(chain.nextvec)-1):
            tmp1 = chain.nextvec[vecindex] 
            tmp2 = chain.nextvec[vecindex+1]
            if (tmp1 != tmp2):
                chain.nextvec[vecindex] = tmp2 
                chain.nextvec[vecindex + 1] = tmp1
            else: 
                # default: do a rigid rotation
                self._do_rigid_rot(vecindex, direction, chain)

        # if possible, 1/3 of the time do a crankshft
        # (1st and 3rd dirs must be different)
        elif (t < 0.66666) and (vecindex < len(chain.nextvec)-2):
            tmp1 = chain.nextvec[vecindex] 
            tmp2 = chain.nextvec[vecindex+2]
            if (t < 0.66666) and (tmp1 != tmp2):
                ### crankshaft move
                chain.nextvec[vecindex] = tmp2 
                chain.nextvec[vecindex + 2] = tmp1
            else: 
                # default: do a rigid rotation
                self._do_rigid_rot(vecindex, direction, chain)

        else: 
            # default: do a rigid rotation
            self._do_rigid_rot(vecindex, direction, chain)

        chain.nextcoords = \
            chain.vec2coords(chain.nextvec, chain.nextcoords)
        chain.nextviable = \
            chain.viability(chain.nextcoords)

    def move3(self, chain):
        """
        Apply moveset 'MC3' to the chain.
        This is just a simple set to change the direction of
        a single chain link.
        Example:
            [0,0,0,0,0] --> [0,0,1,0,0]
        where {0,1,2,3}={n,e,s,w} direction
        About 5% viable moves are expected.
        """
    
        r = self.random.random()
        s = self.random.random()
        vecindex = int(floor((chain.n - 1.0001)*r))

        if s < 0.5:
            direction = 1
        else:
            direction = -1      
  
        chain.nextvec[vecindex] = \
            (chain.nextvec[vecindex] + direction) % 4
        chain.nextcoords = \
            chain.vec2coords(chain.nextvec, chain.nextcoords)
        chain.nextviable = \
            chain.viability(chain.nextcoords)

    def move4(self, chain):
        """
        Apply moveset 'MC4' to the chain:
        This is another vert simple moveset, to just change one angle
        in a rigid rotation.
        Like 'MS3', this generates about 5% viable moves.
        """
        r = random()
        s = random()
        vecindex = int(floor((chain.n - 1.0001)*r))

        if s < 0.5:
            direction = 1
        else:
            direction = -1

        ### a rigid rotation
        self._do_rigid_rot(vecindex, direction, replica)

        chain.nextcoords = \
            chain.vec2coords(chain.nextvec, chain.nextcoords)
        chain.nextviable = \
            chain.viability(chain.nextcoords)

    def metropolis(self, replica):
        """
        Accept Chain.nextvec over Chain.vec according to a Metropolis
        criterion.
        """
        randnum = random()

        # accept with Metroplis criterion
        thisenergy = self.energy(replica.chain) + \
                        self.restraint.energy(replica.chain)  
        boltzfactor = exp((thisenergy - self.lastenergy) / \
                        (self.k*self.temp))

        if randnum < boltzfactor:
            # update the chain
            self._update_chain(replica)
            # update the lastenergy
            self.lastenergy = thisenergy
            return 1
        else:
            return 0

    def _update_chain(self, replica):
        replica.chain.vec[:] = replica.chain.nextvec[:]
        replica.chain.coords[:,:] = replica.chain.nextcoords[:,:]
        replica.chain.viable = replica.chain.nextviable

    def energy(self, chain):
        # """Calculate potential energy of the chain."""
        # num_contacts = 0.0
        # for i, this_H_idx in enumerate(self.H_inds):
        #     for other_H_idx in self.H_inds[i+1:]:
        #         if (other_H_idx - this_H_idx) >= 3:
        #             this_dist = (abs(chain.coords[this_H_idx][0]-chain.coords[other_H_idx][0]) + \
        #                          abs(chain.coords[this_H_idx][1]-chain.coords[other_H_idx][1]))
        #             if this_dist == 1:
        #                 num_contacts += 1
        # return num_contacts * self.epsilon
        return compute_energy(self.epsilon, chain.coords, self.H_inds)


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
        for i in range(0, len(self.contacts)):
            # print 'contact', i, self.contacts[i], chain.coords
            c = self.contacts[i][0]
            d = self.contacts[i][1]
            D = D + (chain.coords[c][0]-chain.coords[d][0])**2
            D = D + (chain.coords[c][1]-chain.coords[d][1])**2
        return D
