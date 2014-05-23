#! /usr/bin/python

import sys
sys.path.append('../')
from os.path import join
from random import Random
from string import strip
from math import exp

from timer import Timer
from hplattice.Config import Config
from hplattice.Monty import randseed
from hplattice.Replica import Replica
from hplattice.Trajectory import Trajectory


CLIST_PATH = '/Users/grollins/work_code/HP_lattice/sequences/clist'
HP_STRING_SET = \
    ((4, ['HPPH', 'HHPH', 'HPHH', 'HHHH'], join(CLIST_PATH, 'hp04')),
     (6, ['HHPPHH', 'HPPHPH', 'HHPHPH', 'HPPHHH'], join(CLIST_PATH, 'hp06')),
     (10, ['HHPPHPPHPH', 'HPHPPHPPHH', 'PHPPHHPPHP', 'HPPHPPHPPH'], join(CLIST_PATH, 'hp10')),
     (12, ['HPHPPHPHPHPH', 'PHPPHPPHPPHP', 'HHHHHPHHPHPH'], join(CLIST_PATH, 'hp12')),
     (14, ['HHHPHHPHHHHPPH', 'HPPPPHPPHPPHPH', 'HHPPHHPHPHHHHH'], join(CLIST_PATH, 'hp14')),
     # (20, ['PHHHPPHHHPPPPPHHPPHP', 'HHHHPPHHHHPHHPHPPHHH',
     # 'HHHPPPPHPPHPPPPHPPHP', 'HHHHPPHHPHHHHHPPHPHH'])
)

g = Random(randseed)
VERBOSE = False
    
##########################
# Functions 

def attemptswap(config, replicas, swaps, viableswaps): 

    # Attempt a swap between replicas    
    if config.SWAPMETHOD == 'random pair':
        # pick pair at random
        r = g.random()
        i = min(int(r * config.NREPLICAS), (config.NREPLICAS - 1))
        j = i
        while j == i:
            s = g.random()
            j = min(int(s * config.NREPLICAS), (config.NREPLICAS - 1))
        # make sure j > i (needed for the swap criterion below)
        if j < i:
            tmp = i
            i = j
            j = tmp
    
    elif config.SWAPMETHOD == 'neighbors':    
        # pick neighboring pair at random 
        r = g.random()
        i = min(int(r * (config.NREPLICAS - 1)), (config.NREPLICAS - 2))
        j = i + 1
                
    else:
        print 'Swap method', config.SWAPMETHOD, 'unknown.  Exiting...'
        sys.exit(1)

    if (VERBOSE):
      print 'REX: Attempting swap between replicas', i, 'and', j
     
    randnum = g.random()

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
        
        retval = 1
        
    else:
        retval = 0

    viableswaps[i] = viableswaps[i] + retval
    viableswaps[j] = viableswaps[j] + retval
    swaps[i] = swaps[i] + 1
    swaps[j] = swaps[j] + 1

    return retval

#####################################
# Main Program
   
def mc_sampling(hp_string, initial_vec, clist_dir):
    # load in config file
    configfile = 'mcrex.conf'
    config = Config( filename=configfile )
    config.HPSTRING = hp_string
    config.INITIALVEC = initial_vec
    config.RESTRAINED_STATE = []
    config.NATIVEDIR = clist_dir

    # look up the native contact state from a pre-compiled flat file     
    if config.STOPATNATIVE == 1:
        nativeclistfile = config.NATIVEDIR + '/' + config.HPSTRING + '.clist'
        fnative = open(nativeclistfile,'r')
        nativeclist = eval(fnative.readline())
        fnative.close()
        if VERBOSE:
            print 'NATIVE CLIST:',nativeclist
    
    
    # Make a list of config.NREPLICAS Replica objects
    replicas = []           # a list of Replica objects
    for i in range(0, config.NREPLICAS):
        replicas.append( Replica(config,i) )
    
    # Each replica is just a container for the objects:
    #     replica[i].chain   <-- the HP chain
    #     replica[i].mc      <-- A Monte Carlo move set object that performs operations on the chain
    
    # Keep track of statistics for the replica exchange simulation, for each replica:

    ### Monte Carlo stats
    steps = []                 # number of total move steps attempted (viable or not)
    viablesteps = []           # number of viable steps
    acceptedsteps = []         # number of viable steps that were accepted under Metropolis criterion 
    moveacceptance = []        # fraction of 
    acceptance = []            # fraction of MC moves accepted

    ### Replica exchange stats 
    swaps = []                 # the number of attemped replica swaps
    viableswaps = []           # the number of viable swaps
    swap_acceptance = []       # the fraction of swaps that were accepted
    
    for i in range(0, config.NREPLICAS):
         steps.append(0)
         viablesteps.append(0)    
         acceptedsteps.append(0)    
         moveacceptance.append(0)    
         acceptance.append(0)    
         swaps.append(0)
         viableswaps.append(0)
         swap_acceptance.append(0)
    
    prodstep = 1
    foundnative = 0
    while (prodstep < (config.MCSTEPS + 1)) & (foundnative == 0):
    
        # Run the replicas for a production cycle...
        for rep in range(0, config.NREPLICAS):

            ### Propose a new MC move
            if strip(config.MOVESET) == 'MS1':
                replicas[rep].mc.move1(replicas[rep])
            elif strip(config.MOVESET) == 'MS2':
                replicas[rep].mc.move2(replicas[rep])
            elif strip(config.MOVESET) == 'MS3':
                replicas[rep].mc.move3(replicas[rep])
            elif strip(config.MOVESET) == 'MS4':
                replicas[rep].mc.move4(replicas[rep])
            else: 
                print 'MC MOVESET=', config.MOVESET, 'not supported!'
                print 'Exiting....'
                sys.exit(1)

            # count this move only if the chain is viable
            if replicas[rep].chain.nextviable == 1:
                ### accept with metroplis criterion
                accept = replicas[rep].mc.metropolis(replicas[rep])
                # if accept is yes, the chain is updated (see Monty.py)
                if (accept):
                    acceptedsteps[rep] = acceptedsteps[rep] + 1
                # keep track of MC steps
                viablesteps[rep] = viablesteps[rep] + 1

            # keep track of total steps
            steps[rep] = steps[rep] + 1

            # HAVE WE FOUND THE NATIVE YET? If so, stop
            if config.STOPATNATIVE == 1:
                thisclist  = replicas[rep].chain.contactstate()
                if (nativeclist == thisclist):
                    foundnative = 1
                    print "Found native at step %d, rep %d" % \
                        (steps[rep], rep)

        # After the production cycle,      
        if (prodstep % config.SWAPEVERY) == 0:

            ### ...after every production run, attempt a SWAP
            success = attemptswap(config, replicas,swaps,viableswaps)

            if (VERBOSE):
              if success:   
                outstring = str(prodstep)+'\tSwap successful!\treplica temps: ['
                for rep in range(0,config.NREPLICAS):
                  outstring = outstring + str(replicas[rep].mc.temp) + ' '
                print outstring +']'

              else:  
                print str(prodstep)+'\tUnsuccessful swap attempt.'

        # Print status
        if (prodstep % config.PRINTEVERY) == 0:

            # calc MC move acceptance
            for rep in range(0,len(steps)):
                if steps[rep] > 0:
                    moveacceptance[rep] = float(viablesteps[rep])/float(steps[rep])
                    if viablesteps[rep] > 0:
                        acceptance[rep] = float(acceptedsteps[rep])/float(viablesteps[rep])
                    else:
                        acceptance[rep] = 0.0
                else:
                    moveacceptance[rep] = 0.0
                    acceptance[rep] = 0.0

            # calc replica swap acceptance
            for rep in range(0,len(steps)):
                if swaps[rep] > 0:
                    swap_acceptance[rep] = float(viableswaps[rep])/float(swaps[rep])
                else:
                    swap_acceptance[rep] = 0.0

            # Output the status of the simulation
            '''
            print prodstep, 'production steps'
            print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s ' % \
                  ('replica','viablesteps','steps','MCaccept','viableswaps',
                   'swaps','SWAPaccept')
            for rep in range(0, len(steps)):
                print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s ' % \
                    (rep,viablesteps[rep],steps[rep],
                     '%1.3f' % moveacceptance[rep], swaps[rep],
                     viableswaps[rep], '%1.3f' % swap_acceptance[rep])
            if config.STOPATNATIVE == 1:
                print 'NATIVE CLIST:', nativeclist
            print '%-8s %-12s %-12s' % \
                ('replica', 'foundnative', 'contact state') 
            for replica in replicas:
                print '%-8d %-12d %s' % \
                    (replica.repnum, foundnative,
                     repr(replica.chain.contactstate()))
            '''

        # Continue to the next production cycle!
        prodstep = prodstep + 1

    # write the acceptance ratios to the <rundir>/data directory
    # faccept = open(config.DATADIR + '/acceptance', 'w')
    # for i in range(0, len(config.REPLICATEMPS)):
    #     tmp = str(config.REPLICATEMPS[i]) + '\t'
    #     tmp = tmp + str(acceptedsteps[i]) + '\t'
    #     tmp = tmp + str(viablesteps[i]) + '\t'
    #     tmp = tmp + str(steps[i]) + '\n'

    #Output the status of the simulation one last time...
    '''
    print prodstep, 'production steps'
    print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s ' % \
            ('replica','viablesteps','steps','MOVEaccept',
             'viableswaps','swaps','SWAPaccept')
    for rep in range(0, len(steps)):
        print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s ' % \
            (rep,viablesteps[rep], steps[rep],
             '%1.3f' % moveacceptance[rep],
              swaps[rep], viableswaps[rep],
              '%1.3f' % swap_acceptance[rep])
    if config.STOPATNATIVE == 1:
        print 'NATIVE CLIST:', nativeclist
    print '%-8s %-12s %-12s' % \
        ('replica','foundnative','contact state')    
    for replica in replicas:
        print '%-8d %-12d %s' % \
            (replica.repnum,foundnative,repr(replica.chain.contactstate()))
    '''

    # faccept.write(tmp)  
    # faccept.close()

    print "last step was", prodstep

def main():
    timing_data = []
    for N, hp_string_list, clist_dir in HP_STRING_SET:
        for this_hp_string in hp_string_list:
            print N, this_hp_string, clist_dir
            initial_vec = [0] * (N - 1)
            with Timer() as t:
                mc_sampling(this_hp_string, initial_vec, clist_dir)
            time_elapsed = t.interval
            this_data = (N, time_elapsed, this_hp_string)
            timing_data.append(this_data)

    for td in timing_data:
        print "%d  %.2e  %s" % td

if __name__ == '__main__':
    main()