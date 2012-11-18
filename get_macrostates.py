#!/usr/bin/env sage -python

import sys
sys.path.append("/home/kenlocey/modules/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
from sage.all import * 
import resource
from os import path, access, R_OK  # W_OK for write permission
import re
sys.path.append("/home/kenlocey/modules/combine_macrostates")
import combine_macrostates as combine

""" These functions analyze data from macroecological and metagenomic datasets. 
    Functions look for text files nameofdataset-data.txt in folder nameofdataset.
    Functions expect the data to be in text file format with the following as
    columns: site, species, abundance of species. For each name in datasets,
    there should be a directory with that name and a text file holding site,
    species, and species abundance data in columns:
	    site1, K sonoriense, 15
	    site1, K baurii, 10
	    site1, K flavescens, 8 
	    site2, K sonoriense, 18
	    site2, K baurii, 13
	    site2, K flavescens, 12 """

""" A list call 'datasets' contains the name of folders, in which, the name of the folder is also a prefix for txt files
that the functions below look for (e.g. folder name is DATA, DATA contains a file DATA-data.txt). """

datasets = ['BBS','CBC','FIA','GENTRY','MCDB','NABC','TERA','AQUA','FUNGI']
#for name in os.listdir('/home/kenlocey/data1'): 
#    datasets.append(name) 


""" Below are some functions to generate figures in recent feasible set manuscript (Locey and White)"""

""" Find random macrostates with the same combinations of N and S as found in datasets.
    The user should check the 'get_macrostates' function in the module 'feasible_functions' to ensure that the 
    desired number of random macrostates, along with other variables, is set. """ 
#ff.get_random_macrostates_for_datasets(datasets)

""" Generate text files of observed vs expected SAD data """
#ff.generate_obs_pred_data(datasets)

""" Generate figures from Locey and White (201?) """
#Figure 1
ff.get_all_SADs([[40,[5,20,35]],[400,[20,90,150]],[1000,[40,140,210]]])

# Figure 2
#ff.obs_pred_r2_multi(datasets)
#ff.plot_obs_pred_sad(datasets)

# Figure 3
#ff.pairwise_r2_obs_feasible(datasets)

# Figure 4
#NS_combos = [[22,10],[22,12],[29,10],[29,12],[32,10],[32,12],[39,10],[39,12],[42,10],[42,12],[49,10],[49,11],[52,10],[53,11],[56,10],[59,10]]
#ff.dataset_NS_combos(NS_combos,'FIA')



""" Generate additional figures and analyses for examining feasible sets and empirical SAD data """
#ff.kdens_feasible_vs_obs(50,10,'FIA')

#ff.kdens_full_feasibles(50,10)

#ff.get_500_RADs([[50,10],[50,20],[60,20],[60,30]])

#ff.plot_obs_exp_evenness(datasets)

#ff.plot_percentile(datasets)

#ff.get_common_combos1(datasets)

""" Combine macrostates that were randomly generated from multiple sources and stored in text files """
#datasets = [] # a list of file names corresponding to folder names that designate datasets
#combine.combine_macrostates(datasets)