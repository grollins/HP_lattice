from string import zfill
from random import random
from math import exp

from .Chain import Chain
from .Monty import Monty


class Replica:
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self, config, repnum, nativeclist=None):
        """Initialize the Replica() object."""
        temp = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.nativeclist = nativeclist
        self.repname = 'rep' + zfill( str(repnum), 2 )
        self.repdir = config.DATADIR + self.repname
        self.chain = Chain(config)
        self.mc = Monty(config, temp, self.chain)
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
        thisclist = self.chain.contactstate()
        return (self.nativeclist == thisclist)


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

    delfactor = ((1/(replicas[j].mc.k*replicas[j].mc.temp)) - (1/(replicas[j].mc.k*replicas[i].mc.temp)))
    delfactor = delfactor * (replicas[j].mc.energy(replicas[j].chain) - replicas[i].mc.energy(replicas[i].chain) )

    boltzfactor = exp(delfactor)
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
