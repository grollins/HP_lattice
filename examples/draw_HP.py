#! /usr/bin/env python

import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import sys
import pickle


fi = open(sys.argv[1],'r')
x = pickle.load(fi)

s = x.pop(1)
print s
s0,s1 = s
HPseq = x.pop(0)
seq = numpy.array(list(HPseq))
H = seq == 'H'
P = seq == 'P'

spring = [s0,s1]
print spring[0]
print spring[1]

for i,conf in enumerate(x):
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    X = []
    Y = []
    for (x,y) in conf:
        X.append(x)
        Y.append(y)
    
    X = numpy.array(X)
    Y = numpy.array(Y)


    ax1.text(X[0]+0.1,Y[0],0,size=10,color='black')
    for j in range(1,len(X)):
        ax1.plot([X[j-1],X[j]],[Y[j-1],Y[j]],'k-')
        ax1.text(X[j]+0.1,Y[j]-0.2,j,size=10,color='black')

    ax1.plot(X[H],Y[H],'bo',markersize=8,markeredgecolor='b')
    ax1.plot(X[P],Y[P],'ro',markersize=8,markeredgecolor='r')

    ax1.set_xlim(min(X)-1,max(X)+1)
    ax1.set_ylim(min(Y)-1,max(Y)+1)
    ax1.grid(True)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    #plot spring
    ax1.plot([X[spring[0]],X[spring[1]]],[Y[spring[0]],Y[spring[1]]],'g--')
    fig.savefig('fig%i.png' % i,format='png')
        

