#! /usr/bin/python

usage = """
enumerate.py <configfile>

Try:  enumerate.py enumerate.conf

This program will read in an HP chain specified in the configure file,
and perform a full enumeration of conformational space.

The problem tablulates:

    1) the density of states (in energies/contacts)
    
    2) the number density of unique contact states, i.e. disjoint collections
       of microscopic conformations all sharing a unique set of interresidue contacts. 

These values are printed as output.

"""


import sys

from Config import *
from Chain import *
from Monty import *
from Replica import *
from Trajectory import *

import numpy
import random
import string
import math
import os

def harmonic(dist,spring,d=1):
    """DO a harmonic potential after d=d, before have no penalty in energy"""
    if dist <= d:
        return 0
    
    return 0.5*spring*(dist-1)


g = random.Random(randseed)


if len(sys.argv) < 2:
    print 'Usage:  enumerate.py <configfile>'
    sys.exit(1)
 
 
VERBOSE = 1
    

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print usage
        sys.exit(1)
        
    configfile = sys.argv[1]
    config = Config( filename=configfile)
    if len(config.RESTRAINED_STATE) > 0:
        restraint = DistRestraint(config.RESTRAINED_STATE,config.KSPRING)
    if VERBOSE: config.print_config()
    
    # create a single Replica
    replicas = [ Replica(config,0) ]
    
    traj = Trajectory(replicas,config)        # a trajectory object to write out trajectories

    nconfs = 0
    contact_states = {}                # dictionary of {repr{contact state}: number of conformations}
    contacts = {}               # dictionary of {number of contacts: number of conformations}
    springs = {}                #dictionary of {distance of restraint: number of restraints}


    #################
    #
    # This is a useful subroutine for enumerating all conformations of an HP chain
    #
    # NOTE: in order for this to work correctly, the initial starting vector must be [0,0,0,....,0]
    # 
    done = 0
    conformations = []
    conformations.append(replicas[0].chain.hpstring)
    if len(config.RESTRAINED_STATE) > 0:
        conformations.append(config.RESTRAINED_STATE[0])
    else:
        conformations.append((0,0))
    while not(done):
                        
        if len(replicas[0].chain.vec) == replicas[0].chain.n-1:    
            if replicas[0].chain.viable:                
                if replicas[0].chain.nonsym():
                    
                    # tally the number of contacts
                    state = replicas[0].chain.contactstate()
                    try:
                        D = restraint.D(replicas[0].chain)
                    except:
                        D = 0

                    ncontacts = len(state)
                    if contacts.has_key((ncontacts,D)) == False:
                        contacts[(ncontacts,D)] = 1
                    else:
                        contacts[(ncontacts,D)] += 1
                    if ncontacts == 7 and D<2:
                       conformations.append(replicas[0].chain.vec2coords(replicas[0].chain.vec))

                    # tally the contact state
                    this_state_repr = repr(state)
                    if contact_states.has_key((this_state_repr,D)) == False:
                        contact_states[(this_state_repr,D)] = 1
                    else:
                        contact_states[(this_state_repr,D)] +=1

                    # tally the number of conformations
                    nconfs = nconfs + 1

                    # write to trajectory
                    if (nconfs % config.TRJEVERY) == 0:
                        traj.queue_trj(replicas[0])
                    # print progress
                    if (nconfs % config.PRINTEVERY) == 0:
                        print '%-4d confs  %s'%(nconfs,replicas[0].chain.vec)
    
                done = replicas[0].chain.shift()
                    
            else:
                done = replicas[0].chain.shift()

        else:
            if replicas[0].chain.viable:
                replicas[0].chain.grow()
            else:
                done = replicas[0].chain.shift()

        if replicas[0].chain.vec[0] == 1:    # skip the other symmetries
            break        
    #
    #
    #################
        
    
    # write the last of the trj and ene buffers
    # and close all the open trajectory file handles
    traj.cleanup(replicas)
    
    # print out the density of contact states
    print
    print 'DENSITY of CONTACT STATES:'
    print '%-40s %s'%('contact state','number of conformations')
    for state in contact_states.keys():
            a,b = state
            #if "(0, 13)" in a:
            #    continue
            #if "(1, 8)" in a:
            #    continue
            #else:
            #if a.count(',') == 9:
            #        print '%-40s %10s %d'%(a,b,contact_states[state])
            print '%-40s %10s %d'%(a,b,contact_states[state])

    # print out the density of states (energies)
    print 
    print 'DENSITY of STATES (in energies/contacts):'
    print '%-20s %-20s %s'%('number of contacts','energy (kT)','number of conformations')
    #Calculate partition function
    Z = 0
    maximum = 0
    fo = open('partition.dat','w')
    N_cont = []
    Rest_dist = []
    Tot_E = []
    Degeneracy = []
    for c in contacts.keys():
        a,b = c
        ener = config.eps*a+harmonic(b,config.KSPRING,d=1.0)
        if a > maximum:
            maximum = a 
            e_max = math.exp(-ener)
        Z += contacts[c]*math.exp(-ener)
        print '%-20d %-20d %-20.2f %d'%(a,b,ener,contacts[c])
        N_cont.append(a)
        Rest_dist.append(b)
        Tot_E.append(ener)
        Degeneracy.append(contacts[c])
    all_vec = [N_cont,Rest_dist,Tot_E,Degeneracy]
    pickle.dump(numpy.array(all_vec),fo)
    fo.close()
    print
    print 'at T = %4.1f K'%config.T
    print e_max,Z,e_max/Z,-0.6*math.log(e_max/Z)

    fo = open('data_conf.dat','w')
    pickle.dump(conformations,fo)
    fo.close
        

    
