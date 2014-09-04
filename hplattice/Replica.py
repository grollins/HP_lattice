from random import random
from math import exp


class Replica(object):
    """
    Each *Replica* is a container for a chain and a monte carlo sampler.
    During replica exchange simulations, the replicas attempt to swap
    temperatures at regular intervals. The success of a swap depends on
    the energies and temperatures of the replicas.

    :param lattice_factory: factory object that knows how to create chains and
                            monty samplers.
    :type lattice_factory: :class:`hplattice.LatticeFactory`
    :param str config: path to configuration file
    :param int repnum: replica number
    :param list nativeclist: optional, native contacts as a list of tuples,
                             example ``[(0, 4), (1, 6)]``
    """
    def __init__(self, lattice_factory, config, repnum, nativeclist=None):
        T = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.nativeclist = nativeclist
        self.chain = \
            lattice_factory.make_chain(config.HPSTRING, config.INITIALVEC)
        self.mc = lattice_factory.make_monty(config, T, self.chain)
        self.mc_move_fcn = self._select_move(config.MOVESET.strip())

    def init_mc_stats(self):
        """
        Initialize stats about viability and acceptance of monte carlo moves.
        """
        # number of total move steps attempted (viable or not)
        self.steps = 0
        # number of viable steps
        self.viablesteps = 0
        # number of viable steps that were accepted under Metropolis criterion 
        self.acceptedsteps = 0
        # fraction of moves that are viable
        self.move_viability = 0
        # fraction of viable moves that are accepted
        self.acceptance = 0

    def record_stats(self, move_is_viable, move_is_accepted):
        """
        Update mc move stats based on outcome of most recent move.

        :param bool move_is_viable: ``True`` if the move was a viable,
                                    non-overlapping chain conformation.
        :param bool move_is_accepted: ``True`` if the move was accepted as the
                                      new conformation of the chain.
        """
        if move_is_viable:
            self.viablesteps += 1
            if move_is_accepted:
                self.acceptedsteps += 1
        self.steps += 1

    def compute_mc_acceptance(self):
        """
        Compute the fraction of mc moves that have been viable and the fraction
        of viable moves that have been accepted.
        """
        if self.steps > 0:
            self.move_viability = \
                (1. * self.viablesteps) / self.steps
            if self.viablesteps > 0:
                self.acceptance = float(self.acceptedsteps)/float(self.viablesteps)
            else:
                self.acceptance = 0.0
        else:
            self.move_viability = 0.0
            self.acceptance = 0.0

    def _select_move(self, move_name):
        ### maps move names to move methods
        if move_name == 'MS1':
            move_fcn = self.mc.move1
        elif move_name == 'MS2':
            move_fcn = self.mc.move2
        elif move_name == 'MS3':
            move_fcn = self.mc.move3
        else: 
            print 'MC MOVESET=', move_name, 'not supported!'
            move_fcn = None
        return move_fcn

    def propose_move(self):
        """
        Do a monte carlo move to produce a new conformation of the chain.

        :return: ``True`` if new conformation is viable.
        :rtype: bool
        """
        self.mc_move_fcn(self.chain)
        return self.chain.nextviable()

    def metropolis_accept_move(self):
        """
        Apply metropolis criterion to determine whether new conformation of
        chain should replace current conformation.

        :return: ``True`` if new conformation should be accepted.
        :rtype: bool
        """
        return self.mc.metropolis(self)

    def is_native(self):
        """
        Check if current contacts match contacts of native state.

        :return: ``True`` if current contacts match native contacts.
        :rtype: bool
        """
        this_contact_list = self.contactstate()
        return (self.nativeclist == this_contact_list)

    def contactstate(self):
        """
        Get contacts of current chain conformation.

        :return: list of tuples, example ``[(0, 4), (1, 6)]``
        :rtype: list
        """
        return self.chain.contactstate()

    def get_vec(self):
        """
        Get chain vectors.

        :return: chain vectors
        :rtype: :class:`hplattice.Chain.Chain.Vectors`
        """
        return self.chain.vec

    def get_T(self):
        """
        Get current temperature.

        :return: temperature
        :rtype: float
        """
        return self.mc.temp

    def energy(self):
        """
        Compute energy of current chain conformation.

        :return: energy
        :rtype: float
        """
        return self.mc.energy(self.chain)

    def kT(self):
        """
        :return: :math:`k_b * T`
        :rtype: float
        """
        return self.mc.kT()


def attemptswap(swap_method, replicas):
    """
    Attempt swap of replicas.

    :param str swap_method: ``'random pair'`` to randomly choose two replicas
                            to swap; ``'neighbors'`` to randomly choose one
                            replica ``i`` and swap it with its ``i+1`` neighbor.
    :param list replicas: list of :class:`Replica` objects
    :return: the indices of the two replicas and the success or failure of the swap.
    :rtype: (int, int, bool)
    """
    N = len(replicas)
    # Attempt a swap between replicas
    if swap_method == 'random pair':
        # pick pair at random
        r = random()
        i = min(int(r * N), (N - 1))
        j = i
        while j == i:
            s = random()
            j = min(int(s * N), (N - 1))
        # make sure j > i (needed for the swap criterion below)
        if j < i:
            tmp = i
            i = j
            j = tmp
    
    elif swap_method == 'neighbors':    
        # pick neighboring pair at random 
        r = random()
        i = min(int(r * (N - 1)), (N - 2))
        j = i + 1
                
    else:
        print 'Swap method', swap_method, 'unknown.'
 
    randnum = random()

    ### if proposing i-->j, 
    ### 
    ###   del = d(beta_j - beta_i) * d(E_j - E_i)
    ###
    ### (see Hansmann 1997)

    boltzfactor = _compute_boltz_factor(replicas[i], replicas[j])

    if randnum < boltzfactor:

        # swap the ***temperatures***
        temp_i = replicas[i].mc.temp
        temp_j = replicas[j].mc.temp
        tempfromrep_i = replicas[i].mc.tempfromrep
        tempfromrep_j = replicas[j].mc.tempfromrep
        replicas[i].mc.temp = temp_j
        replicas[j].mc.temp = temp_i
        replicas[i].mc.tempfromrep = tempfromrep_j
        replicas[j].mc.tempfromrep = tempfromrep_i
        
        swap_success = True
        
    else:
        swap_success = False

    return i, j, swap_success

def _compute_boltz_factor(replica_i, replica_j):
    ### compute boltzmann factor for a pair of replicas
    delfactor = (1. / replica_j.kT()) - (1. / replica_i.kT())
    delfactor = delfactor * (replica_j.energy() - replica_i.energy())
    boltzfactor = exp(delfactor)
    return boltzfactor
