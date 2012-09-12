#!/usr/bin/env sage -python

import sys
sys.path.append("/home/kenlocey/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
from sage.all import * 
import resource
from os import path, access, R_OK  # W_OK for write permission
import re

""" These functions analyze data from macroecological and metagenomic datasets. They compare observed SADs
    to that expected from sampling the feasible set and generate figures of the Locey and White (2013) paper.
    Functions look for text files nameofdataset-data.txt in folder nameofdataset. Functions expect the data
    to be in text file format with the following as columns: site, species, abundance of species.
    For each name in datasets, there should be a directory with that name and a text file holding site,
	species, and species abundance data in columns:
	    site1, K sonoriense, 15
	    site1, K baurii, 10
	    site1, K flavescens, 8 
	    site2, K sonoriense, 18
	    site2, K baurii, 13
	    site2, K flavescens, 12 
    Uncomment function calls and print statements to execute."""

#datasets = ['BBS','CBC','FIA','GENTRY','MCDB','NABC','TERA','AQUA','FUNGI']
#datasets = ['COAL_PROD','COAL_CONS','COAL_EMIT','OIL_PROD','OIL_CONS','OIL_EMIT','NGAS_PROD','NGAS_CONS','NGAS_EMIT']
datasets = ['WHEAT_PROD','WHEAT_SUP','WHEAT_WAST','RICE_PROD','RICE_SUP','RICE_WAST','CORN_PROD','CORN_SUP','CORN_WAST']
#datasets = ['OIL_PROD','OIL_CONS','OIL_EMIT']
ff.get_macrostates(datasets) 
""" The above function finds random macrostates with the same combinations of N and S as found in datasets.
    The user should check the get_macrostates function in the module feasible_functions to ensure that the 
    desired number of random macrostates, along with other variables, is chosen. """ 

# Figure 1
#ff.get_all_SADs(50,10)
#print 'figure 1 complete'

# Figure 2
#ff.Evar_kdens_feasible_vs_obs(50,10,'FIA')
#print 'figure 2 complete'

# Figure 3
#ff.Evar_kdens_full_feasibles(50,10)
#print 'figure 3 complete'

#Figure 4
#ff.get_500_RADs([[50,10],[50,20],[60,20],[60,30]])
#print 'figure 4 complete'

# Figure 5
#ff.generate_obs_pred_data(datasets)
#ff.obs_pred_r2_multi(datasets)
#ff.plot_obs_pred_sad(datasets)
#print 'figure 5 complete'

#Figure 6
#ff.pairwise_r2_obs_feasible(datasets)
#print 'figure 6 complete'

# Figure 7
#ff.plot_obs_exp_evenness(datasets)
#print 'figure 7 complete'

# Figure 8
#ff.plot_percentile(datasets)
#print 'figure 8 complete'

#Figure 9
#datasets = ['BBS','CBC','FIA','GENTRY','MCDB','NABC','CATLIN','CHU','LAUB','HYDRO']
#ff.get_common_combos(datasets)
#NS_combos = [[22,10],[22,12],[29,10],[29,12],[32,10],[32,12],[39,10],[39,12],[42,10],[42,12],[49,10],[49,11],[52,10],[53,11],[56,10],[59,10]]
#ff.FIA_NS_combos(NS_combos,'FIA')
#print 'figure 9 complete'




