#!/usr/bin/python

import sys
sys.path.append("/your_path/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
import resource
from os import path, access, R_OK  # W_OK for write permission
import re

""" These functions analyze data from macroecological and metagenomic datasets.
    They reproduce figures and analyses from Locey and White (2013).
    The user should check that paths, options, and variables in the functions
    called from feasible_functions. Uncomment to run functions. """

datasets = ['BBS','CBC','FIA','GENTRY','MCDB','TERA','AQUA','FUNGI']

""" Find random macrostates for combinations of N and S found in datasets. 
    The user should check the 'get_random_macrostates_for_datasets' function in the module
    'feasible_functions' to ensure that the desired number of random macrostates, along
    with other options are chosen.
    
    Alternatively, the user can email Ken Locey (ken@weecology.org) to request a copy
    of the macrostates we generated. Copies will be sent ASAP. """ 
#ff.get_random_macrostates_for_datasets(datasets)

""" Generate text files of observed vs expected SAD data. Copies of the files used in
    Locey and White (2013) can be found in the 'public_data' folder of this repository """
#ff.generate_obs_pred_data(datasets)

""" Generate figures from Locey and White (201?) """
# generate figure 1
ff.get_all_SADs([[1000,[40,140,210]]])

# generate figure 2
ff.mode_evenness_skew()

# generate figure 3 and get overall r2 values for each dataset
ff.obs_pred_r2_multi(datasets)
ff.plot_obs_pred_sad(datasets)

# generate figure 4
ff.pairwise_r2_obs_feasible(datasets)

#ff.get_common_combos2(['FIA'],20,0,0,20)
""" This function is not necessary to reproduce results from Locey and White (2013), as the list
    of N-S combinations used in that paper are provided below. However, if the user wishes to generate figure 5
    for a different set of N-S combos, they should run the above function to find which combos meet the criteria.
    
    Arguments for the above function:
    1. Name of dataset. This corresponds to a file ending in '_data.txt', e.g. 'FIA' matches FIA_data.txt
    2. Minimum N. This is the smallest value of N to be considered.
    3. Largest difference between values of N; allows for similar but not identical values of N. When set at
       zero, all values of N must be identical
    4. Largest difference between values of S; allows for similar but not identical values of S. When set at 
       zero, all values of S must be identical
    5. Minimum sample size. N-S combos with less than this number of replicates are ignored,
       with or without a tolerance for varying N and S.

    Below: NS_combos was derived using the above function. """    
NS_combos = [[22,10],[22,12],[29,10],[29,12],[32,10],[32,12],[39,10],[39,12],[42,10],[42,12],[49,10],[49,11],[52,10],[53,11],[56,10],[59,10]]
# generate figure 5
ff.dataset_NS_combos(NS_combos,'FIA')
