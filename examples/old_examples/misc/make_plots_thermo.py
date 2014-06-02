#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import thermodynamics
import math
import numpy


def plot_E_vs_degen(S):
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)

    for sys in S.keys():
        temp = S[sys]
        logg = numpy.array([math.log(g) for g in temp['g']])
        E =numpy.array(temp['E'])
        ind = E.argsort()
        ax1.plot(logg[ind],E[ind])

    fig.savefig('plot.png',format='png')

def plot_properties(S,prop):
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)

    x = numpy.arange(len(S.keys()))
    y = [] 
    labels = []
    for sys in S.keys():
        temp = S[sys]
        y.append(temp[prop])
        labels.append(sys)

    ax1.bar(x,y,color='w',edgecolor='k')
    ax1.set_xticks(x+0.5)
    ax1.set_xticklabels(labels)
    print labels

    fig.savefig('%s.png' % prop,format='png')
 






canonical = thermodynamics.thermo('no_springs/partition.dat')
good_spring = thermodynamics.thermo('good_ECO/partition.dat')
bad_spring = thermodynamics.thermo('bad_ECO/partition.dat')
best_spring = thermodynamics.thermo('good_spring/partition.dat')
worst_spring = thermodynamics.thermo('bad_spring/partition.dat')

systems = {}
systems['canonical'] = canonical
systems['good_spring'] = good_spring
systems['bad_spring'] = bad_spring
systems['best_spring'] = best_spring
systems['bad2_spring'] = worst_spring

plot_E_vs_degen(systems)
plot_properties(systems,'Z')
plot_properties(systems,'U')
plot_properties(systems,'F')
plot_properties(systems,'S1')
plot_properties(systems,'Cv')
