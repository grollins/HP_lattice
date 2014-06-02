#! /usr/bin/env python

import numpy
import sys
import pickle
import math

KT = 0.6
T = 300.

def partition_function(E,degeneracy):
    '''calculate the partition function'''
    z = numpy.array([d*math.exp(-ener) for d,ener in zip(degeneracy,E)])
    return z.sum()

def averages(param,E,deg,z):
    '''
        calculate the average of a parameter: 
        <E> = sum(E*deg*exp{-E/KT}/z
    '''
    av = numpy.array([p*d*math.exp(-m) for p,d,m in zip(param,deg,E)])
    return av.sum()/z

def entro(E,deg,z):
    p = numpy.array([d*math.exp(-m) for d,m in zip(deg,E)])
    p = p/z
    s = numpy.array([pr*math.log(pr) for pr in p])
    return -KT/T*s.sum()

def correct_degeneracy(e,g):
    '''
    When we add springs, the energies make it so that states with different amount
    of contacts have the same energy. I was keeping track of these as separate, have to unite them if
    we are going to plot density of states,...
    '''

    track = {}
    for ene,d in zip(e,g):
        if track.has_key(ene) == False:
            track[ene] = d
        else:
            track[ene] += d

    e = []
    g = []
    for key in track.keys():
        e.append(float(key))
        g.append(int(track[key]))
    return numpy.array(e),numpy.array(g)

def thermo(filen):

    fi = open(filen,'r')
    data = pickle.load(fi)

    contacts = data[0,:]
    restraint= data[1,:]
    energy = data[2,:]
    degeneracy = data[3,:]

    energy,degeneracy = correct_degeneracy(energy,degeneracy)

    z = partition_function(energy,degeneracy)

    E_av = averages(energy*KT,energy,degeneracy,z)

    E_sq_av = averages(energy*energy*KT*KT,energy,degeneracy,z)

    Cv_heat_capacity = (E_sq_av-E_av*E_av) / (KT*T)

    A_free_energy = -KT*math.log(z)

    S_entropy = (E_av - A_free_energy)/T

    S_method_2 = entro(energy,degeneracy,z)

    obj = {}
    obj['E'] = energy
    obj['g'] = degeneracy
    obj['Z'] = z
    obj['F'] = A_free_energy
    obj['U'] = E_av
    obj['S1'] = S_entropy
    obj['S2'] = S_method_2
    obj['Cv'] = Cv_heat_capacity

    print obj['Z'],obj['F'],obj['U'],obj['S1'],obj['S2'],obj['Cv']

    return obj


if __name__ == "__main__":
    system = thermo(sys.argv[1])

    print system['Z'],system['F'],system['U'],system['S1'],system['S2'],system['Cv']

