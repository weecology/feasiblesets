#!/usr/bin/sage

""" Code for integer partitioning generated by Justin Kitzes. """

'''
Routines for calculating central tendency of all possible integer partitions of 
n (ie, N or individuals) with k (ie, S for species) parts.
'''

from __future__ import division
from sage.all import *
import numpy as np
import matplotlib.pyplot as plt


def psi_dist(abund, N, S, set=None):
    '''
    Calculate the probability of finding 0 to S species with given abundance.

    Parameters
    ----------
    abund : int
        Abundance at which to calculate distribution, from 0 to N (technically 
        can only be >0 up to N - S - 1).
    N : int
        Number of total individuals in community
    S : int
        Number of species in community
    set : int
        Optional set of allowable parts in partition

    Returns
    -------
    psi : ndarray
        Array giving probability of i species, from 0 to S, with abundance 
        abund.
    '''

    tol = 1e-6  # Tolerance for ending loop

    if set == None:
        parts = number_of_partitions(N, S)
    else:
        parts = number_of_partitions_restricted(N, set, S)

    psi_cum = np.zeros(S+2)

    for i in xrange(0, S+2):

        if N - abund*i < 0 or S - i < 0:
            break
        else:
            if set == None:
                psi_cum[i] = number_of_partitions(N-abund*i, S-i) / parts
            else:
                psi_cum[i] = number_of_partitions_restricted(N-abund*i, set, 
                                                             k=S-i) / parts

        if psi_cum[i] < tol:  # If low probability of higher abunds
            break

    return psi_cum[:-1] - psi_cum[1:]


def samp_part(N, S):
    '''
    Get uniform random sample of integer partitions of N with length S.

    Parameters
    ----------
    N : int
        Number of total individuals in community
    S : int
        Number of species in community

    Returns
    -------
    part : ndarray
        Random partition of length S.
    '''

    N_run = N
    S_run = S

    part = np.zeros(S)
    pind = 0

    for abund in xrange(1, N+1):  # Start at lowest abund = 1

        # Calculate psi for this abund, with only addends >= abund allowed
        psi = psi_dist(abund, N_run, S_run, set=range(abund, N+1))
        cum_psi = np.cumsum(psi)

        # Choose num spp w/ this abund
        num_choice = np.argmax(cum_psi > np.random.random())
        part[pind:(pind+num_choice)] = abund

        # Increment counters
        pind += num_choice
        N_run -= num_choice * abund
        S_run -= num_choice

        if S_run == 0:  # If have chosen abund for all species
            break
            
    part = part.tolist()    
    part.reverse()
    SAD = []
    for p in part:SAD.append(int(p))
    return SAD


#path_save = '/home/kenlocey/'
#a = psi_dist(2,1000,100)
#plt.close()
#plt.plot(a)
#plt.savefig(path_save + 'JKitzes_fig.png')
