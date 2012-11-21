#! /usr/bin/env python

import os

state=['no_springs','good_spring','bad_spring','good_ECO','bad_ECO']
spring=['RESTRAINED_STATE        []','RESTRAINED_STATE       [(0,11)]','RESTRAINED_STATE        [(0,7)]','RESTRAINED_STATE        [(6,13)]','RESTRAINED_STATE        [(1,8)]']
temp=['200','250','300','350','400','450','500']
epsilon=['eps   -1.5','eps  -1.2','eps  -1.0','eps  -0.8571','eps  -0.75','eps  -0.6667','eps  -0.6']
k=['KSPRING 1.5','KSPRING 1.2','KSPRING 1.0','KSPRING 0.8571','KSPRING 0.75','KSPRING 0.6667','KSPRING 0.6']


for i in range(len(temp)):
    if not os.path.exists(temp[i]):
        os.mkdir(temp[i])
    os.chdir(temp[i])

    for j in range(len(state)):
        if not os.path.exists(state[j]):
            os.mkdir(state[j])
        os.chdir(state[j])
        fo = open('enumerate_long.conf','w')
        string ='''HPSTRING                HHPHHPHHHPPHHH
        INITIALVEC              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        randseed                345
        NREPLICAS               1
        REPLICATEMPS            [%s]
        EXPDIR                  ./enumerate_data
        PRINTEVERY              1000
        TRJEVERY                1000
        ENEEVERY                1000
        %s
        %s
        %s''' % (temp[i],spring[j],epsilon[i],k[i])

        fo.write(string)
        fo.close() 
        os.system('../../enumerate_restraints.py enumerate_long.conf > enumerate_long.log')
        os.system('../../draw_HP.py data_conf.dat')
        os.system('../../thermodynamics.py partition.dat')
        os.chdir('..')

    os.chdir('..')

