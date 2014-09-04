from random import random
from math import floor, exp
from numpy import array, int32


BOLTZ_CONST = 0.001987  # (kcal/K.mol) Boltzmann's constant


class Monty(object):
    """
    A collection of functions to perform Monte Carlo move-set operations
    on an HP chain.

    :param config: configuration parameters for chain and simulation
    :type config: :class:`hplattice.Config.Config`
    :param float temp: temperature (K)
    :param chain: do monte carlo on this chain
    :type chain: :class:`hplattice.Chain.Chain`
    """

    def __init__(self, config, temp, chain):
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

        # (both of these are copies from Config() )
        self.restraint = DistRestraint(config.RESTRAINED_STATE, config.KSPRING)
        self.lastenergy = self.energy(chain) + self.restraint.energy(chain)

    def kT(self):
        """
        :return: :math:`k_b * T`
        :rtype: float
        """
        return BOLTZ_CONST * self.temp

    def move1(self, chain, vecindex=None, direction=None):
        """
        Apply moveset MC1 (Dill and Chan, 1994, 1996) to the chain:
        
        1. three-bead flips
        2. rigid rotations

        :param chain: apply move to this chain
        :type chain: :class:`hplattice.Chain.Chain`
        :param int vecindex: optional, vector to move. will be chosen randomly
                             if no value is specified.
        :param int direction: optional, ``1`` for clockwise, ``-1`` for
                              counterclockwise. will be chosen randomly if no
                              value is specified.
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
        Apply moveset MC2 (Dill and Chan, 1994, 1996) to the chain:
        
        1. three-bead flips
        2. crankshaft moves
        3. rigid rotations

        :param chain: apply move to this chain
        :type chain: :class:`hplattice.Chain.Chain`
        :param int vecindex: optional, vector to move. will be chosen randomly
                             if no value is specified.
        :param int direction: optional, ``1`` for clockwise, ``-1`` for
                              counterclockwise. will be chosen randomly if no
                              value is specified.
        :param float moveseed: optional, *moveseed* :math:`<1/3` for three-bead
                             flip;  :math:`1/3<=` *moveseed* :math:`<2/3` for
                             crankshaft; :math:`2/3<=` *moveseed* for rigid
                             rotation.
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
        Apply rigid rotations to the chain:

        1. rigid rotations

        :param chain: apply move to this chain
        :type chain: :class:`hplattice.Chain.Chain`
        :param int vecindex: optional, vector to move. will be chosen randomly
                             if no value is specified.
        :param int direction: optional, ``1`` for clockwise, ``-1`` for
                              counterclockwise. will be chosen randomly if no
                              value is specified.
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
        Judge the next conformation of the chain according to Metropolis
        criterion: :math:`e^{-\Delta E/kT}`.

        :param replica: The replica containing the chain that should be judged.
        :type replica: :class:`hplattice.Replica.Replica`
        :return: ``True`` if next conformation should be accepted.
        :rtype: bool
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
            return True
        else:
            return False

    def energy(self, chain):
        E, state = chain.energy(self.epsilon)
        return E


class DistRestraint:
    """
    A harmonic constraint based on squared distance :math:`D = d^2`,
    where :math:`D = \sum_{i,j} d^2_{ij}` over all contacts.

    :param list contacts: list of tuples, example ``[(0, 4), (1, 6)]``
    :param float kspring: spring constant for restraint
    """

    def __init__(self, contacts, kspring):
        self.contacts = contacts
        self.kspring = kspring # (J/[lattice space]^2)

    def energy(self, chain):
        """
        Compute the restraint energy of a given chain.

        :param chain: Compute the restraint energy of this chain.
        :type chain: :class:`hplattice.Chain.Chain`
        :return: energy of the distance restraint
        :rtype: float
        """
        return self.kspring * self.D(chain)

    def D(self, chain):
        """
        Compute the sum of squared-distances of a given chain.

        :param chain: Compute the sum of squared-distances of this chain.
        :type chain: :class:`hplattice.Chain.Chain`
        :return: the sum of squared-distances over the selected contacts.
        :rtype: float
        """
        D = 0.0
        coords = chain.get_coord_array()
        for i in range(0, len(self.contacts)):
            # print 'contact', i, self.contacts[i], chain.coords
            c = self.contacts[i][0]
            d = self.contacts[i][1]
            D = D + (coords[c][0]-coords[d][0])**2
            D = D + (coords[c][1]-coords[d][1])**2
        return D
