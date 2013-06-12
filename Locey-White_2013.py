#!/usr/bin/python

import sys
sys.path.append("/home/kenlocey/modules/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
import resource
from os import path, access, R_OK  # W_OK for write permission
import re
sys.path.append("/home/kenlocey/modules/combine_macrostates")
import combine_macrostates as combine

""" These functions analyze data from macroecological and metagenomic datasets.
    They reproduce figures and analyses from Locey and White 2013.
    The user should check the function called from feasible_functions and 
    combine_macrostates for correct file paths, variable values, ect. """

#datasets = []
#for name in os.listdir('/home/kenlocey/data'): 
#    datasets.append(name) 
datasets = ['BBS','CBC','FIA','GENTRY','MCDB','TERA','AQUA','FUNGI']

""" Find random macrostates with the same combinations of N and S as found in datasets.
    The user should check the get_random_macrostates_for_datasets function in the module
    feasible_functions to ensure that the desired number of random macrostates, along
    with variable values is chosen. """ 
#ff.get_random_macrostates_for_datasets(datasets)

""" Generate text files of observed vs expected SAD data """
#ff.generate_obs_pred_data(datasets)

""" Generate figures from Locey and White (201?) """
# Figure 1
ff.get_all_SADs([[1000,[40,140,210]]])

# Figure 2
#ff.mode_evenness_skew()

# Figure 3
#ff.obs_pred_r2_multi(datasets)
#ff.plot_obs_pred_sad(datasets)

# Figure 4
#ff.pairwise_r2_obs_feasible(datasets)

# Figure 5
#ff.get_common_combos2(['FIA'],20,0,0,20)
#NS_combos = [[22,10],[22,12],[29,10],[29,12],[32,10],[32,12],[39,10],[39,12],[42,10],[42,12],[49,10],[49,11],[52,10],[53,11],[56,10],[59,10]]
#ff.dataset_NS_combos(NS_combos,'FIA')
