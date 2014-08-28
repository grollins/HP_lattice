from string import zfill
from random import random
from math import exp


class Replica(object):
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self, lattice_factory, config, repnum, nativeclist=None):
        """Initialize the Replica() object."""
        T = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.nativeclist = nativeclist
        self.repname = 'rep' + zfill( str(repnum), 2 )
        self.repdir = config.DATADIR + self.repname
        self.chain = \
            lattice_factory.make_chain(config.HPSTRING, config.INITIALVEC)
        self.mc = lattice_factory.make_monty(config, T, self.chain)
        self.mc_move_fcn = self._select_move(config.MOVESET.strip())

    def init_mc_stats(self):
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
        if move_is_viable:
            self.viablesteps += 1
            if move_is_accepted:
                self.acceptedsteps += 1
        self.steps += 1

    def compute_mc_acceptance(self):
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
        if move_name == 'MS1':
            move_fcn = self.mc.move1
        elif move_name == 'MS2':
            move_fcn = self.mc.move2
        elif move_name == 'MS3':
            move_fcn = self.mc.move3
        elif move_name == 'MS4':
            move_fcn = self.mc.move4
        else: 
            print 'MC MOVESET=', move_name, 'not supported!'
            move_fcn = None
        return move_fcn

    def propose_move(self):
        self.mc_move_fcn(self.chain)
        return self.chain.nextviable()

    def metropolis_accept_move(self):
        return self.mc.metropolis(self)

    def is_native(self):
        this_contact_list = self.contactstate()
        return (self.nativeclist == this_contact_list)

    def contactstate(self):
        return self.chain.contactstate()

    def get_vec(self):
        return self.chain.vec

    def get_T(self):
        return self.mc.temp

    def energy(self):
        return self.mc.energy(self.chain)

    def kT(self):
        return self.mc.kT()

def attemptswap(swap_method, replicas):
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
    delfactor = (1. / replica_j.kT()) - (1. / replica_i.kT())
    delfactor = delfactor * (replica_j.energy() - replica_i.energy())
    boltzfactor = exp(delfactor)
    return boltzfactor
