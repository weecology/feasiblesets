#!/usr/bin/env sage -python

from sage.all import *
import csv
import sys
sys.path.append("/home/kenlocey/modules/pymods")
import macroecotools
import os
from os import path, access, R_OK  # W_OK for write permission
import  matplotlib.pyplot as plt
from pylab import *
import numpy as np
from scipy.stats import gaussian_kde
from scipy import stats
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import random
from random import choice
import re
from multiprocessing import Pool, freeze_support
import itertools
from operator import itemgetter
from decimal import *
import math
import random, decimal

########################################################################################################
######   A Section devoted to evenness indices and descriptive statistical functions ###################

def Shannons_H(sad):
    H = 0
    for i in sad:
        p = float(i)/float(sum(sad))
        H += p*ln(RDF(p))
    return H*-1
    
def simplest_gini(x):
    """Return computed Gini coefficient of inequality. This function was found at http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py

	:note: follows basic formula
	:see: `calc_gini2`
	:contact: aisaac AT american.edu
	"""
	x = sorted(x)  # increasing order
	n = len(x)
	G = sum(xi * (i+1) for i,xi in enumerate(x))
	G = 2.0*G/(n*sum(x)) #2*B
	return G - 1 - (1./n)
    
def gini_sample(SADs):
    """ Compute Gini's coefficient for each macrostate in a random sample """
    Gs = []
    for sad in SADs:
        G = simplest_gini(sad)
        Gs.append(G)
    return Gs


def Mcintosh_evenness(SAD):
    S = len(SAD)
    N = sum(SAD)
    sum_n = 0
    for n in SAD: sum_n += n**2
    U = np.sqrt(sum_n)    
    E = (N - U)/(N - (N/np.sqrt(S)))
    return E
def Mcintosh_sample(SADs):
    """ Compute McIntosh evenness for each macrostate in a random sample """
    Es = []
    for sad in SADs:
        E = Mcintosh_evenness(sad)
        Es.append(E)
    return Es

def pielous_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0
    for p in SAD:
        H += -(p/N)*ln(RDF(p/N)) 
    J = H/ln(RDF(S))
    return J
def pielous_sample(SADs):
    """ Compute Pielous evenness for each macrostate in a random sample """
    Js = []
    for sad in SADs:
        J = pielous_evenness(sad)
        Js.append(J)
    return Js    


def NHC_evenness(SAD):
    x_list = range(1,len(SAD)+1)
    y_list = np.log(SAD)
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)
    return slope
def NHC_sample(SADs):
    """ Compute NHC evenness for each macrostate in a random sample """
    NHCs = []
    for sad in SADs:
        NHC = NHC_evenness(sad)
        NHCs.append(NHC)
    return NHCs


def Heips_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0
    for p in SAD:
        H += -(p/N)*ln(RDF(p/N)) 
    H = (np.exp(H) - 1)/(S - 1)
    return H
def Heips_sample(SADs):
    """ Compute Heips evenness for each macrostate in a random sample """
    Hs = []
    for sad in SADs:
        H = Heips_evenness(sad)
        Hs.append(H)
    return Hs


def simpsons_evenness(SAD):
    D = 0.0
    N = float(sum(SAD))
    S = len(SAD)
    for x in SAD:
        D += (x*(x-1))/(N*(N-1))
    E = (1/D)/S
    return E
def simpsons_sample(SADs):
    """ Compute Simpsons evenness for each macrostate in a random sample """
    Es = []
    for sad in SADs:
        E = simpsons_evenness(sad)
        Es.append(E)
    return Es


def berger_parker(SAD):
    bp = float(max(SAD))/sum(SAD)
    return bp
def berger_parker_sample(SADs):
    """ Compute Berger-Parkers evenness index for each macrostate in a random sample """
    BPs = []
    for sad in SADs:
        BP = berger_parker(sad)
        BPs.append(BP)
    return BPs

    
def EQ_evenness(SAD):
    S = len(SAD)
    y_list = np.log(SAD)
    x_list = []
    for rank in range(1,S+1): x_list.append(rank/float(S))
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)
    Eq = (-2.0/pi)*arctan(slope)
    return Eq
def EQ_sample(SADs):
    """ Compute EQ evenness index for each macrostate in a random sample """
    EQs = []
    for sad in SADs:
        Eq = EQ_evenness(sad)
        EQs.append(Eq)
    return EQs   


def e_var(SAD):
    P = np.log(SAD)
    S = len(SAD)
    X = 0
    for x in P:
        X += (x - mean(P))**2/S
    evar = 1 - 2/pi*arctan(X) 
    return(evar)
def Evars_sample(unique_sads):
    """ Compute Evar for each macrostate in a random sample """
    Evars = []
    for sad in unique_sads:
        Evar = e_var(sad)
        Evars.append(Evar)
    return Evars


def get_skews(_list):

    skews = []
    for i in _list:
        skews.append(stats.skew(i))
    
    return skews


def get_modes(_list,which):

    modes = []
    for i in _list:
        _mode = mode(i)
        if which is 'high':
            modes.append(int(max(_mode)))            
        if which is 'low':
            modes.append(int(min(_mode)))
    
    return modes


def get_modal(_list):
    
    """ Finds the mode from a kernel density function across a sample """
    exp_mode = 0.0
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    density.covariance_factor = lambda : .001
    density._compute_covariance()
    D = [xs,density(xs)]
    d = 0
    maxd = 0.0
    while d < len(D[1]):
        if D[1][d] > maxd:
            maxd = D[1][d]
            exp_mode = D[0][d]
        d += 1
    return exp_mode
    
def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
    
def get_kdens(_list):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : 0.5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

 
#######################################################################################################
#### A section devoted to finding constraint combinations and empirical DOWs/SADs from datasets #######

def get_expSADs_fromfile(dataset):
    
    PATH = '/home/kenlocey/data/' + dataset + '/' + dataset + '_obs_pred.txt'
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
        DATA = open(PATH,'r')
        ct1 = 0
        ct2 = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        SADs = []
        
        for d in DATA:
            ct1+=1
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 1 and sum(SAD) <= 100000:
                    SAD.sort()
                    SAD.reverse()
                    SADs.append(SAD) # can also append, site_name, len(SAD), and sum(SAD)
                    ct2+=1
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
        DATA.close()
        return(SADs)
    

def get_SADs(dataset):

    DATA = open('/home/kenlocey/data/' + dataset + '/' + dataset + '-data.txt','r')
    ct1 = 0
    ct2 = 0
    d = DATA.readline()
    m0 = re.match(r'\A\S*',d).group()
    m2 = int(re.findall(r'\d*\S$',d)[0])
    SAD = [int(m2)]
    SADs = []
        
    for d in DATA:
        ct1+=1
        m1 = re.match(r'\A\S*',d).group()
        if m1 == m0:
            m2 = int(re.findall(r'\d*\S$',d)[0])
            if m2 > 0:
                SAD.append(m2)
                
        else:
            site_name = m0
            m0 = m1
            if len(SAD) > 9 and sum(SAD) <= 100000:
                SAD.sort()
                SAD.reverse()
                SADs.append(SAD) # can also append, site_name, len(SAD), and sum(SAD)  !!THIS NEEDS TO BE DEALT WITH!! !!THIS IS GOING TO TRIP SOMEBODY UP!!
                ct2+=1
            SAD = []
            abundance = int(re.findall(r'\d*\S$',d)[0])
            if abundance > 0:SAD.append(abundance)
    DATA.close()
    return(SADs)

def get_NS_combos_labels(datasets):
    
    NS_combos = []
    NS_combos_dataset = []
    for dataset in datasets:
        DATA = open('/home/kenlocey/data/' + dataset + '/' + dataset + '-data.txt','r')
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        for d in DATA:
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 9:
                    NS_combos.append([sum(SAD),len(SAD)])
                    NS_combos_dataset.append([sum(SAD),len(SAD),dataset])
                    
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
        DATA.close()
        
    unique_NS_combos = [list(x) for x in set(tuple(x) for x in NS_combos)]
    return [unique_NS_combos,NS_combos_dataset]
    
    
def get_NS_combos(datasets):
    NS_combos = []
    total_combos = 0
    for dataset in datasets:
        SADs = get_SADs(dataset)
        print len(SADs),'usable sites in',dataset      
        for SAD in SADs:
            NS_combos.append([sum(SAD),len(SAD)])
        
    NS_combos = [list(x) for x in set(tuple(x) for x in NS_combos)]
    print len(NS_combos),'unique NS_combos' 
    
    return (NS_combos)


def get_all_partitions(N,S):
    
    parts = []
    for p in Partitions(N,length=S):
        parts.append(list(p))

    return parts    
        
        
########################################################################################################
######   A Section devoted to finding macrostates/integer partitions ############################

def random_parts(N,S,size): # A newly discovered method for generating random samples of feasible sets
    
    SADs = []
    while len(SADs) < size:
        SAD = list(Partitions(N).random_element())
        
        while len(SAD) != S:
            if len(SAD) == 1:break
            SAD = list(Partition(SAD).conjugate())
            if len(SAD) == S or len(SAD) < 3:break
            r1 = choice(list(set(SAD)))
            SAD.remove(r1)            
            r2 = choice(SAD)
            SAD.remove(r2)
            SAD = list(Partition(SAD).conjugate()) # get the conjugate before appending
            SAD.append(r1+r2)
            SAD.sort()
            SAD.reverse()
            if len(SAD) < 3:break
                
        if len(SAD)==S:
            SADs.append(SAD)
            #print len(SADs),N,S
            
    SADs = [list(x) for x in set(tuple(x) for x in SADs)]
    return SADs
    
def worker2(NS_combo):
    """thread worker function"""
    set_random_seed()
    random_macros = random_parts(NS_combo[0],NS_combo[1],63)
    return random_macros    

def get_random_macrostates(NS_combo):    
    N = int(NS_combo[0])
    S = int(NS_combo[1])
    ct = 0
    rand_macros = []
    while len(rand_macros) < 63:
        macro = Partitions(N).random_element()
        if len(macro) == S:
            rand_macros.append(macro)
        ct+=1                
    rand_macros = [list(x) for x in set(tuple(x) for x in rand_macros)]
    return rand_macros

def worker1(NS_combo):
    """thread worker function"""
    set_random_seed()
    random_macros = get_random_macrostates(NS_combo)
    return random_macros
    
def get_rand_sample(NS_combo): #choose between worker2 (random partitioning alg derived by KJL) and worker1 (random partitioning alg provided by SAGE)
    
    unique_SADs = []
    pool = Pool()
    unique_SADs = pool.map(worker2, [NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo])
    """ worker1 and worker2 call different functions for generating random macrostates. worker1 uses Sage's function (def. worker2 uses the function developed by Ken Locey (faster)."""
    pool.close()
    pool.join()
    return unique_SADs


def get_random_macrostates_for_NScombos(NS_combos):    
    while NS_combos:
        ct = 0
        for NS_combo in NS_combos:
            print len(NS_combos),'NS combinations left'
            ct+=1
            N = int(NS_combo[0])
            S = int(NS_combo[1])
            OUT = open('/home/kenlocey/combined1/' + str(N) + '-' + str(S) + '.txt','a+')
            macros = len(set(OUT.readlines()))
            p = float(number_of_partitions(N,S))/float(number_of_partitions(N))
            if p > 10.0**-6 and macros < 400 and macros != number_of_partitions(N,S):
                rand_macros = get_rand_sample(NS_combo) # Use multiprocessing
                for i in rand_macros:
                    for p in i:
                        print>>OUT,p
                OUT.close()
                NS_combos.remove(NS_combo)
                    
            elif len(NS_combos) == 1:
                NS_combos.remove(NS_combo)
                OUT.close
                break
            else:    
                NS_combos.remove(NS_combo)
                OUT.close()
            if not NS_combos:break
    return

def get_random_macrostates_for_datasets(datasets):

    NS_combos = get_NS_combos(datasets)
    get_random_macrostates_for_NScombos(NS_combos)
    
    return

""" below are two flavors for getting common combinations of N and S from datasets.
    The first only simply prints them to the screen. The second all for greater
    specification and returns them in list form."""
   
def get_common_combos1(datasets):
    n = 0
    all_NS_combos = get_NS_combos(datasets)
    unique_NS_combos = [list(x) for x in set(tuple(x) for x in all_NS_combos)]
    common_combos = []
    for combo in unique_NS_combos:
        if all_NS_combos.count(combo) >= S and combo[0] >= N:
            num = all_NS_combos.count(combo)
            print combo,num
            n+=num            
            common_combos.append(combo)
    print 'number of unique N-S combos:',len(common_combos),'  number of sites:',n


def get_common_combos2(datasets,Nmin,Ndiff,Sdiff):
    
    combos = get_NS_combos_labels(datasets)
    unique_NS_combos = combos[0]
    NS_combos_dataset = combos[1]
    ct3 = 0
    common_combos = []
    for combo1 in unique_NS_combos:
        ct = 0
        _list = []
        for combo2 in NS_combos_dataset:
            if combo1[0] >= Nmin:
                if abs(combo2[0] - combo1[0]) <= Ndiff:
                    if abs(combo2[1] - combo1[1]) <= Sdiff:
                        ct+=1
                        _list.append(combo2[2])
        
        u_list = list(set(_list))
        if len(u_list) >= 5:
            ct3+=1
            #print combo1,' ',ct,' ',len(u_list)
            common_combos.append([combo1[0],combo1[1],ct,u_list])                        
            #if ct3 == 3:break    
    return (common_combos)
    


########################################################################################################
##### A Section devoted to examining randoms samples of feasibles sets and empirical data ##############

    
def get_hottest_SAD(unique_SADs):
    """ Find the SAD in a random sample with the greatest average commonness 
        among its ranked abundance states. This SAD is taken to represent the 
        central tendency of the set, based on the SAD shape. """
    
    if len(unique_SADs) > 500:
        unique_SADs = random.sample(unique_SADs,500)
    
    N = sum(unique_SADs[0])
    S = len(unique_SADs[0])
    a1 = 0 # SAD mean
    v1 = 0 # SAD variance
    for rad in unique_SADs:
        in_common = []
        ct1 = 0
        for a in rad: # for each rank
            c = 0
            for sad in unique_SADs: 
                if a == sad[ct1]:
                    c += 1
            in_common.append(ln(RDF(c)))
            ct1 += 1
        a2 = mean(in_common)
        v2 = variance(in_common)  
        if a2 > a1:
            a1 = a2
            v1 = v2
            xRAD = rad
        elif a2 == a1:
            if v2 < v1:
                a1 = a2
                v1 = v2
                xRAD = rad
    #percentile_evar = stats.percentileofscore(sample_evar,obs_evar)
    return xRAD
                
def hottest_SAD_full_feasible(N,S):
    """ Find the SAD in the feasible set with the greatest average commonness 
        among its ranked abundance states. This SAD is taken to represent the 
        central tendency of the feasible set. """ 
    
    unique_SADs = get_all_partitions(N,S)
    expected_SAD = get_hottest_SAD(unique_SADs)

    return expected_SAD
     
def generate_obs_pred_data(datasets):
    
    NS_combos = []
    for dataset in datasets:
        
        """ To prevent the script from tripping on the last line, make the last line of the datafile any sequence of non-data
            related alphanumeric characters (i.e. 0X0)
        """    
        OUT1 = open('/home/kenlocey/data/'+dataset+'/'+ dataset + '_obs_pred.txt','w')
        OUT1.close()
        DATA = open('/home/kenlocey/data/'+dataset+'/'+ dataset + '-data.txt','r')
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        SADs = []
        
        for d in DATA:
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:
                    SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 9 and sum(SAD) < 100000:
                    SADs.append((site_name,sum(SAD),len(SAD),SAD))
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:
                    SAD.append(abundance)
        print dataset,len(SADs)
        num = 0
        for site in SADs:
            site_name = site[0]	
            N = site[1]
            S = site[2]
            SAD = site[3]
            if max(SAD) > 1:
                
                unique_SADs = []
                PATH = '/home/kenlocey/combined1/'+str(N)+'-'+str(S)+'.txt'
                if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                    data = open(PATH,'r')
                    macrostates = data.readlines()
                    for macrostate in macrostates:
                        sad = eval(macrostate)
                        unique_SADs.append(sad)
                    data.close()
                
                if len(unique_SADs) < 300:
                    continue
                if len(unique_SADs) > 500:
                    unique_SADs = random.sample(unique_SADs,500)
                unique_SADs = [list(x) for x in set(tuple(x) for x in unique_SADs)]
                num += 1
                NS_combos.append([N,S])
                #continue
                expSAD = get_hottest_SAD(unique_SADs) # The expected SAD from the random sample  
                SAD.sort()
                SAD.reverse()
                
                #r2 = macroecotools.obs_pred_rsquare(np.log10(SAD), np.log10(expSAD))
                #print dataset,N,S,site_name,r2,' ',berger_parker(SAD),berger_parker(expSAD),' ',simplest_gini(SAD),simplest_gini(expSAD),' ',e_var(SAD),e_var(expSAD) 
                
                ct = 0
                OUT1 = open('/home/kenlocey/data/'+dataset+'/'+ dataset + '_obs_pred.txt','a')
                while ct < len(expSAD): # write to file, by cite, observed and expected ranked abundances
                    print>>OUT1, site_name,SAD[ct],expSAD[ct]
                    ct += 1
                OUT1.close()
        
        DATA.close()
        print dataset,' ',num,'sites\n'  
    num_NScombos = len([list(x) for x in set(tuple(x) for x in NS_combos)])
    print 'NS combos: ',num_NScombos    

def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "S15,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    return data


def hist_mete_r2(sites, obs, pred):  # TAKEN FROM Macroecotools or the mete_sads.py script used for White et al. (2012)
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in sites:
        obs_site = obs[sites==site]
        pred_site = pred[sites==site]
        r2 = macroecotools.obs_pred_rsquare(obs_site, pred_site)
        r2s.append(r2)
    hist_r2 = np.histogram(r2s, range=(0, 1))
    xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
    xvals = xvals[0:len(xvals)-1]
    yvals = hist_r2[0]
    plt.plot(xvals, yvals, 'k-', linewidth=2)
    plt.axis([0, 1, 0, 1.1 * max(yvals)])
    
        

def obs_pred_r2_multi(datasets, data_dir='/home/kenlocey/data/'): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)' 
    for i, dataset in enumerate(datasets):
        obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt') 
        obs = ((obs_pred_data["obs"]))
        pred = ((obs_pred_data["pred"]))
        print dataset,' ',macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))     


def plot_obs_pred_sad(datasets, data_dir='/home/kenlocey/data/', radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    """Multiple obs-predicted plotter""" 
    fig = plt.figure()
    for i, dataset in enumerate(datasets):
        print dataset
        obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt') 
        site = ((obs_pred_data["site"]))
        obs = ((obs_pred_data["obs"]))
        pred = ((obs_pred_data["pred"]))
        
        axis_min = 0.5 * min(obs)
        axis_max = 2 * max(obs)
        ax = fig.add_subplot(3,3,i+1) 
        macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1, 
                                            plot_obj=plt.subplot(3,3,i+1))      
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        #plt.subplots_adjust(left=0.2, bottom=0.12, right=0.8, top=0.92, wspace=0.29, hspace=0.21)  
        plt.subplots_adjust(wspace=0.5, hspace=0.3)
        # Create inset for histogram of site level r^2 values
        axins = inset_axes(ax, width="30%", height="30%", loc=4)
        hist_mete_r2(site, np.log10(obs), np.log10(pred))
        plt.setp(axins, xticks=[], yticks=[])
        
    plt.savefig(dataset + 'obs_pred_plots.png', dpi=800, bbox_inches = 'tight', pad_inches=0)  
       

   
def kdens_full_feasibles(N,S):
    """ Plot kernel density curves of Evar for macrostates of feasible sets
        based on different N and S """
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    parts = get_all_partitions(N,S)
    Evars = Evars_sample(parts) # This could be a sample based on another metric (above)
    D = get_kdens(Evars)
    plt.xlim(0.0, 1.0)
    plt.plot(D[0],D[1],color='black',lw=5)
    
    parts = get_all_partitions(N,S+10)
    Evars = Evars_sample(parts) # This could be a sample based on another metric (above)
    D = get_kdens(Evars)
    plt.xlim(0.0, 1.0)
    plt.plot(D[0],D[1],color='gray',lw=5)
    plt.axvline(x=0.673,ymin=0,ymax=10,color='black',ls='--',lw=3) # plot a vertical line at the mode
    plt.setp(ax, xticks=[0.2,0.6,1.0],yticks=[0,2,4,6])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    
    plt.savefig('Figure3.png', dpi=400, pad_inches=0)     


def kdens_feasible_vs_obs(N,S,dataset):
    """ Plot kernel density curves of Evar for macrostates of a feasible set
        and a sample of macrostates """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    SADs = get_SADs(dataset)
    
    for sad in SADs:
        if len(sad) != S and sum(sad) != N:
            SADs.remove(sad)
    
    print len(SADs),'usable sites in',dataset
    DATA.close()
    
    Evars = Evars_sample(SADs) # This could be a sample based on another metric (above)
    D = get_kdens(Evars)
    D = get_kdens(N,S)
    plt.plot(D[0],D[1],color='black',lw=5)
    
    Evars = Evars_sample(SADs) # This could be a sample based on another metric (above)
    D = get_kdens(Evars)
    D = get_kdens_obs(SADs)
    plt.plot(D[0],D[1],color='gray',lw=5)
    
    plt.setp(ax, xticks=[0.2,0.6,1.0],yticks=[1,3,5])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    plt.savefig('Figure2.png', dpi=400, pad_inches=0) 

   

def histOutline(dataIn, *args, **kwargs): # found at http://www.scipy.org/Cookbook/Matplotlib/UnfilledHistograms
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)
    
    stepSize = binsIn[1] - binsIn[0]
    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]
    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0
    
    return (bins, data)


def get_all_SADs(NSlist): # Figure 1 Locey and White (2013)        ##########################################################################################################
    
    
    i = 1
    fig = plt.figure()
    for _list in NSlist:
        print i
        N = _list[0]
        Slist = _list[1]
        if N > 10:
            for S in Slist:
                get_random_macrostates_for_NScombos([[N,S]])
        
        ct = 0
        ax = fig.add_subplot(3,3,i)
        for S in Slist:
            clr = 0
            if ct == 0: clr= '#1E90FF' #blue 
            elif ct == 1: clr= '0.35'  #grey
            elif ct == 2: clr= '#FF34B3' #red
            
            if N > 10:
                parts = []
                PATH = '/home/kenlocey/combined/'+str(N)+'-'+str(S)+'.txt'
                data = open(PATH,'r')
                macrostates = data.readlines()
                for macrostate in macrostates:
                    part = eval(macrostate)
                    parts.append(part)
                data.close()
            else:
                parts = get_all_partitions(N,S)
            skews = get_skews(parts)
            D = get_kdens_choose_kernel(skews,0.5)
            plt.plot(D[0],D[1],color = clr,lw=3, alpha = 0.99)
            ct+=1
        plt.setp(ax, xticks=[-2,0,2,4,6], yticks=[0.2,0.4,0.6])
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.xlabel("Skewnness",fontsize=10)
        plt.ylabel("pdf",fontsize=10)
        i+=1
        
        
        ct = 0
        ax = fig.add_subplot(3,3,i)
        for S in Slist:
            clr = 0
            if ct == 0: clr= '#1E90FF' #blue 
            elif ct == 1: clr= '0.35'  #grey
            elif ct == 2: clr= '#FF34B3' #red
            macros = []
            if N > 10:
                PATH = '/home/kenlocey/combined/'+str(N)+'-'+str(S)+'.txt'
                data = open(PATH,'r')
                macrostates = data.readlines()
                for macrostate in macrostates:
                    macro = eval(macrostate)
                    macros.append(macro)
                data.close()
            else:
                macros = get_all_partitions(N,S)
            
            modes = []
            for macro in macros:
                _mode = mode(macro)
                modes.extend(_mode)
            x = [0,1,2,3,4,5,6,7,8,9]
            y = [0]*10
            for m in modes:
                if m == 1: y[0]+=1
                elif m == 2: y[1]+=1
                elif m <= 4: y[2]+=1
                elif m <= 8: y[3]+=1
                elif m <= 16: y[4]+=1
                elif m <= 32: y[5]+=1
                elif m <= 64: y[6]+=1
                elif m <= 128: y[7]+=1
                elif m <= 256: y[8]+=1
                elif m <= 512: y[9]+=1
            
            Y = [0]*10
            cty = 0
            for g in y:
                Y[cty] = float(g)/float(sum(y))
                cty += 1
            
            plt.plot(x,Y,'-o',color=clr, lw=2) #label= 'N='+str(N)+' S='+str(S)
            ct+=1
            
        plt.xlim(-1.0,7)
        #plt.ylim(0,0.3)
        #plt.yscale('log')
        plt.xlabel("log2(modal abundance)",fontsize=10)
        plt.ylabel("frequency",fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=7)
        #plt.legend(loc=1,prop={'size':8})
        plt.setp(ax, xticks=[0,1,2,3,4,5,6,7])
        i+=1
        
        
        ct = 0
        ax = fig.add_subplot(3,3,i)
        for S in Slist:
            clr = 0
            if ct == 0: clr= '#1E90FF' #blue 
            elif ct == 1: clr= '0.35'  #grey
            elif ct == 2: clr= '#FF34B3' #red
            parts = []
            PATH = '/home/kenlocey/combined/'+str(N)+'-'+str(S)+'.txt'
            data = open(PATH,'r')
            macrostates = data.readlines()
            for macrostate in macrostates:
                part = eval(macrostate)
                parts.append(part)
            data.close()
            
            if ct == 0 or ct == 2: plt.bar([0],[0], color=clr, linewidth=0, label= 'N='+str(N)+' S='+str(S))
            if ct == 1:
                for part in parts:
                    x = [0,1,2,3,4,5,6,7,8,9]
                    y = [0]*10
                    for p in part:
                        if p >= 512: y[9]+=1
                        elif p >= 256: y[8]+=1
                        elif p >= 128: y[7]+=1
                        elif p >= 64: y[6]+=1
                        elif p >= 32: y[5]+=1
                        elif p >= 16: y[4]+=1
                        elif p >= 8: y[3]+=1
                        elif p >= 4: y[2]+=1
                        elif p >= 2: y[1]+=1
                        elif p == 1: y[0]+=1
                    plt.bar(x,y, color=clr, linewidth=0, align='center', alpha = 0.07)
                plt.bar([0],[0], color=clr, linewidth=0, align='center', label= 'N='+str(N)+' S='+str(S))      
            ct+=1           
            plt.xlim(-1,8)
            plt.xlabel("log2(abundance)",fontsize=10)
            plt.ylabel("frequency",fontsize=10)
            plt.tick_params(axis='both', which='major', labelsize=7)
            plt.setp(ax, xticks=[0,2,4,6,8])
            plt.legend(loc=1,prop={'size':7})
        i+=1
        
    plt.subplots_adjust(wspace=0.35, hspace=0.35)    
    plt.savefig('Fig1-'+str(N)+'-'+str(S)+'.png', dpi=400, pad_inches=0)
    print 'done'
    
    
def plot_obs_exp_evenness(datasets):
    fig = plt.figure()
    i = 0
    for dataset in datasets:
        print dataset
        SADs = []
        OBS_evar = []
        EXP_evar = []
        EXP_EQ = []
        OBS_EQ = []
        OBS_NHC = []
        EXP_NHC = []
        
        SADs = get_SADs(dataset)
        for SAD in SADs:
            N = sum(SAD)
            S = len(SAD)
            PATH = '/home/kenlocey/combined/' + str(N) + '-' + str(S) + '.txt'
            if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                data = open(PATH,'r')
                macros = set(data.readlines())
                if len(macros) >= 100:
                    if len(macros) > 500:
                        macros = random.sample(macros,50)
                    OBS_evar.append(e_var(SAD))
                    OBS_E.append(simpsons_evenness(SAD))

                    OBS_NHC.append(NHC_evenness(SAD))
                    
                    macrostates = []
                    for macro in macros:
                        macrostates.append(eval(macro))
                    expSAD = get_hottest_SAD(macrostates) # The expected SAD from the random sample  
                    if len(expSAD)!=len(SAD) or sum(expSAD)!=sum(SAD):
                        print 'SADs differ in N or S'
                        sys.exit()
                    
                    EXP_evar.append(e_var(expSAD))
                    EXP_E.append(simpsons_evenness(expSAD))

                    EXP_NHC.append(NHC_evenness(expSAD))
                    
                    data.close()
                else:data.close()
            
        print len(OBS_evar),'usable sites'
        ax = fig.add_subplot(3,3,i+1)
        
        r2 = macroecotools.obs_pred_rsquare(np.array(OBS_E),np.array(EXP_E)) # r-squared for the 1:1 line (obs vs. exp)
        slope,intercept,r_value,p_value,std_err = stats.linregress(EXP_E,OBS_E)
        print 'r-squared for obs-v-pred Simpsons evenness:',r2,r_value**2
        plt.scatter(EXP_E, OBS_E, c='b',marker='o',lw=0,s=10,alpha= 0.3)
        
        r2 = macroecotools.obs_pred_rsquare(np.array(OBS_evar),np.array(EXP_evar)) # r-squared for the 1:1 line (obs vs. exp)
        slope,intercept,r_value,p_value,std_err = stats.linregress(EXP_evar,OBS_evar)
        print 'r-squared for obs-v-pred Evar:',r2,r_value**2
        plt.scatter(EXP_evar, OBS_evar, c='r',marker='o',lw=0,s=10,alpha=0.4)
              
        r2 = macroecotools.obs_pred_rsquare(np.array(OBS_NHC),np.array(EXP_NHC)) # r-squared for the 1:1 line (obs vs. exp)
        slope,intercept,r_value,p_value,std_err = stats.linregress(EXP_NHC,OBS_NHC)
        print 'r-squared for obs-v-pred NHC evenness:',r2,r_value**2
        
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.plot([0,1],[0,1], 'k-')	
        plt.setp(ax, xticks=[0.2,0.4,0.6,0.8], yticks=[0.2,0.4,0.6,0.8])
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.subplots_adjust(wspace=0.3, hspace=0.3)
        
        # Create inset for kdens of site level comparisons of evenness to the feasible set
        axins = inset_axes(ax, width="40%", height="40%", loc=2)
        axis_min = min(EXP_NHC)
        if min(OBS_NHC) < axis_min: axis_min = min(OBS_NHC)
        axis_max = max(EXP_NHC)
        if max(OBS_NHC) > axis_max: axis_max = max(OBS_NHC)
        
        plt.xlim(axis_min,axis_max)
        plt.ylim(axis_min,axis_max)
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.scatter(EXP_NHC, OBS_NHC, c='0.25',marker='o',lw=0,s=10,alpha=0.3)
        ##ticks = np.arange(axis_min,axis_max, abs(axis_max-axis_min)/3.0)
        ##ticks = np.round_(ticks,decimals=2)
        ##plt.tick_params(axis='both', which='major', labelsize=8)
        plt.setp(axins, xticks=[], yticks=[])
        
        i+=1
    plt.savefig('obsEVEN_expEVEN_plots.png', dpi=400, bbox_inches = 'tight', pad_inches=0.1)  
        
        
def plot_percentile(datasets):
    fig = plt.figure()
    i = 1
    
    ax = fig.add_subplot(1,1,i)
    
    #SADs = []
    #x_list_ginis = []
    #y_list_ginis = []
    
    for dataset in datasets:
        
        x_list_ginis = []
        y_list_ginis = []
        c1 = decimal.Decimal(str(random.random()))
        c2 = decimal.Decimal(str(random.random()))
        c3 = decimal.Decimal(str(random.random()))
        
        SADs = get_SADs(dataset)
        print dataset,len(SADs)
        
        if dataset == 'BBS': dataset = 'Breeding Birds'
        elif dataset == 'CBC': dataset = 'Land Birds'
        elif dataset == 'NABC': dataset = 'Butterflies'
        elif dataset == 'HYDRO': dataset = 'microbe metagenome'
        elif dataset == 'AQUA': dataset = 'microbe metagenome'
        elif dataset == 'TERA': dataset = 'microbe metagenome'
        elif dataset == 'CATLIN': dataset = 'microbe metagenome'
        elif dataset == 'CHU': dataset = 'microbe metagenome'
        elif dataset == 'FUNGI': dataset = 'Fungal metagenome'
        elif dataset == 'MCDB': dataset = 'Mammals'
        elif dataset == 'FIA' or dataset == 'GENTRY': dataset = 'Trees'
        elif dataset == 'LAUB': dataset = 'microbe metagenome'
        elif dataset == 'OIL_PROD':dataset = 'Oil production'
        elif dataset == 'OIL_EMIT':dataset = 'CO2 emissions, Oil'
        elif dataset == 'NGAS_PROD':dataset = 'Natural Gas production'
        elif dataset == 'NGAS_EMIT':dataset = 'CO2 emissions, Nat. Gas'
        elif dataset == 'WHEAT_SUP':dataset = 'Wheat supply'
        elif dataset == 'WHEAT_WAST':dataset = 'Waste of wheat'
        elif dataset == 'RICE_SUP':dataset = 'Rice supply'    
        elif dataset == 'RICE_PROD':dataset = 'Rice production'
        elif dataset == 'rural_pop':dataset = 'Rural population'
        elif dataset == 'urban_pop':dataset = 'Urban population'
        elif dataset == 'ag_pop':dataset = 'Agricultural population'
        elif dataset == 'pelagic_sup':dataset = 'Pelagic fish supply'
        elif dataset == 'carbamate use':dataset = 'Carbamate use'
        elif dataset == 'cereals_prod':dataset = 'Cereal grains produced'
        elif dataset == 'Org_P_pesticide':dataset = 'Organophosphates use'
        elif dataset == 'NGAS_CONS':dataset = 'Nat. Gas consumption'
        elif dataset == 'OIL_CONS':dataset = 'Oil consumption'
        elif dataset == 'land_dev_GCS':dataset = 'Land development, GCS'
        elif dataset == 'land_dev_NCS':dataset = 'Land development'
        elif dataset == 'livestock_assets_NCS':dataset = 'Livestock assets'
        elif dataset == 'livestock_assets_GCS':dataset = 'Livestock assets, GCS'
        elif dataset == 'COAL_CONS':dataset = 'Coal consumption'
        elif dataset == 'CORN_SUP':dataset = 'Corn supply'
        
        if len(SADs) > 100:
            SADs = random.sample(SADs,100)                    
        
        for SAD in SADs:
            N = sum(SAD)
            S = len(SAD)
            PATH = '/home/kenlocey/combined1/' + str(N) + '-' + str(S) + '.txt'
            if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                data = open(PATH,'r')
                macros = set(data.readlines())
                if len(macros) >= 100:
                    if len(macros) > 500:
                        macros = random.sample(macros,500)
                    obs_gini = simplest_gini(SAD)
                    x_list_ginis.append(obs_gini)
                    macrostates = []
                    for macro in macros:
                        macrostates.append(eval(macro))
                    sample_gini = []
                    for macro in macrostates:
                        sample_gini.append(simplest_gini(macro))
                        
                    percentile_gini = stats.percentileofscore(sample_gini,obs_gini)
                    #print len(sample_gini),' ',percentile_gini
                    y_list_ginis.append(percentile_gini)
                    data.close()
                    
                else:data.close()
        
        if len(y_list_ginis) > 0: #and len(SADs) > 30:
            #slope,intercept,r_value,p_value,std_err = stats.linregress(x_list_ginis,y_list_ginis)
            #print 'percentile Gini: r-value:',r_value,'p-value:',p_value,'slope:',slope,' r-squared:',r_value**2
            #m,b = np.polyfit(x_list_ginis,y_list_ginis,1)
            #plt.plot(x_list_ginis, np.array(x_list_ginis)*m +b, c='b',lw=2) 
            plt.scatter(x_list_ginis,y_list_ginis,color=(c1,c2,c3),marker='o',s=100,lw=0,alpha=0.6)
            plt.plot([-10],[-10],color=(c1,c2,c3), label=dataset,lw=4) 
            
    #plt.scatter(x_list_ginis,y_list_ginis,color=(c1,c2,c3),marker='o',s=30,lw=0,alpha=0.5)
        
    plt.xlim(0,1.6)
    plt.ylim(0,100)
    plt.setp(ax, xticks = [0.0,0.2,0.4,0.6,0.8,1.0], yticks=[20,40,60,80,100])
    plt.tick_params(axis='both', which='major', labelsize=14)
    leg = plt.legend(loc=1,prop={'size':14})
    leg.draw_frame(False)
      
    # Create inset for kdens of site level comparisons of evenness to the feasible set
    #slope,intercept,r_value,p_value,std_err = stats.linregress(x_list_NHC,y_list_NHC)
    #print 'percentile NHC: r-value:',r_value,'p-value:',p_value,'slope:',slope,' r-squared:',r_value**2,'\n'
    #axins = inset_axes(ax, width="40%", height="40%", loc=2)
    #plt.scatter(x_list_NHC, y_list_NHC, c='0.25',marker='o',lw=0,s=10,alpha=0.8)
    #plt.setp(axins, xticks=[], yticks=[])
        
    i+=1
    plt.xlabel("Gini's inequality",fontsize=18)
    plt.ylabel("Percentile of a random sample",fontsize=18)    
    #plt.legend(loc=1,prop={'size':10})
    plt.savefig('percentile.png', dpi=1000, bbox_inches = 'tight', pad_inches=0.1)  
    
    
def get_500_RADs(NS_combos): 
    
    i = 1
    fig = plt.figure()
    for combo in NS_combos:
        
        ax = fig.add_subplot(2,2,i)
        N = int(combo[0])
        S = int(combo[1])
        ct = 0
        SADs = []
        while ct < 500:
            macro = Partitions(N).random_element()
            if len(macro) == S:
                ct+=1
                SADs.append(list(macro))
                (bins, n) = histOutline(list(macro))
                plt.plot(bins, n, 'k-',lw=1,alpha=0.05)
                
        exp_SAD = get_hottest_SAD(SADs)        
        (bins, n) = histOutline(exp_SAD)
        plt.plot(bins, n, 'r-',lw=2)
        if S == 10:
            plt.setp(ax,xticks=[1,10,20,30],yticks=[5,10])
            plt.ylim(0,10)
            plt.xlim(0,30)
        if S == 20 and N == 60:
            plt.setp(ax,xticks=[1,5,15,25],yticks=[5,10,15,20])
            plt.ylim(0,20)
            plt.xlim(0,25)
        if S == 20 and N == 50:
            plt.setp(ax,xticks=[1,5,10,15,20],yticks=[5,10,15,20])
            plt.ylim(0,20)
            plt.xlim(0,20)    
        if S == 30:
            plt.setp(ax,xticks=[1,10,20,30],yticks=[10,20,30])
            plt.ylim(0,30)
            plt.xlim(0,25)
        #locs,labels = xticks()
        #xticks(locs, map(lambda x: int(x+1), locs))
        plt.subplots_adjust(left=0.2, bottom=0.12, right=0.8, top=0.92, wspace=0.29, hspace=0.21)  
                
        # Create inset for Evar
        axins = inset_axes(ax, width="50%", height="50%", loc=1)
        expEvar = e_var(exp_SAD)
        plt.axvline(x=expEvar,ymin=0,ymax=10,color='red',ls='--',lw=1.5)
        D = get_kdens_obs(SADs)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='black',lw=3.8)
        D = get_kdens(N,S)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='red',lw=1.3)
        plt.setp(axins, xticks=[0.2,0.8], yticks=[])
        
        i+=1
        print 'finished:'+str(N)+' '+str(S)
        
    plt.savefig('Figure4.png', dpi=400, pad_inches=0)


    
def dataset_NS_combos(NS_combos,dataset):
    """ Plot kernel density curves of species evenness for macrostates of feasible sets
        based on different N and S """
    
    i = 1
    fig = plt.figure()
    for combo in NS_combos:
        ax = fig.add_subplot(4,4,i)
        N = int(combo[0])
        S = int(combo[1])
        DATA = open('/home/kenlocey/data/' + dataset + '/' + dataset + '-data.txt','r')
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        SADs = []
        
        for d in DATA:
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) == S and sum(SAD) == N:
                    SAD.sort()
                    SAD.reverse()
                    SADs.append(SAD)
                    
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
        
        DATA.close()
        _vector = Evars_sample(SADs) # This can also be berger_parker_sample(), Heips_sample(), mode_sample(), etc.
        D = get_kdens(_vector)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='0.6',lw=3) #  blue is '#1E90FF', red '#FF34B3'
        
        rand_macros = []
        _numparts = 0
        while _numparts < 500:
            macro = list(Partitions(N).random_element())
            if len(macro) == S:
                rand_macros.append(macro)
                _numparts+=1
                            
        SADs = [list(x) for x in set(tuple(x) for x in rand_macros)]
        print N,S,len(SADs)
        Evars = Evars_sample(SADs)
        D = get_kdens(Evars)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='k',lw=3) # 
        
        maxd = 0.0
        xspot = 0
        d = 0
        while d < len(D[1]):
            if D[1][d] > maxd:
                maxd = D[1][d]
                xspot = D[0][d]
            d += 1
        plt.axvline(x=xspot,ymin=0,ymax=10,color='black',ls='--',lw=1.5) # plot a vertical line at the mode
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.setp(ax, xticks=[0.2,0.4,0.6,0.8])
        plt.subplots_adjust(wspace=0.35, hspace=0.35)
        i+=1  
    plt.savefig('FIA_feasible_vs_obs.png', dpi=400, pad_inches=0)     
    


def compare(size,datasets):
    
    NS_combos = get_NS_combos(datasets)
    print len(NS_combos)
    ct = 0
    
    for combo in NS_combos:
        
        N = int(combo[0])
        S = int(combo[1])
        
        if N > 20000:continue
        PATH = '/home/kenlocey/combined/'+str(N)+'-'+str(S)+'.txt'
        if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
            fig = plt.figure()
            
            SADs  = []
            data = open(PATH,'r')
            sampleSADs = data.readlines()
            if len(sampleSADs) < 500:continue
            for sad in sampleSADs:
                SADs.append(eval(sad))  
            if len(SADs) > 500:
                SADs = random.sample(SADs,500)
            
            _vector = gini_sample(SADs) # This can also be berger_parker_sample(), Heips_sample(), mode_sample(), etc.
            D = get_kdens_obs(_vector)
            plt.xlim(0.0, 1.0)
            plt.plot(D[0],D[1],color='black',lw=2)    
            
            SADs = []
            SADs = random_parts(N,S,size)
            _vector = gini_sample(SADs) # This can also be berger_parker_sample(), Heips_sample(), mode_sample(), etc.
            D = get_kdens_obs(_vector)
            plt.xlim(0.0, 1.0)
            plt.plot(D[0],D[1],color='red',lw=2)
            data.close()
            
            print 'finished:'+str(N)+' '+str(S)
            plt.savefig('/home/kenlocey/ZPICS/N'+str(N)+'-S'+str(S)+'-'+str(size)+'macros.png', dpi=400, pad_inches=0)
        else: print 'no file for '+str(combo[0])+'-'+str(combo[1])



def pairwise_r2_obs_feasible(datasets):
    
    i = 1
    fig = plt.figure()
    for dataset in datasets:
        ax = fig.add_subplot(3,3,i)
        P_errs = []
        SADs = get_SADs(dataset)
        ct = 0

        for SAD in SADs:
            N = sum(SAD)
            S = len(SAD)
            PATH = '/home/kenlocey/combined1/' + str(N) + '-' + str(S) + '.txt'
            if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                data = open(PATH,'r')
                macros = set(data.readlines())
                if len(macros) >= 400:
                    if len(macros) >= 500:
                        macros = random.sample(macros,500)
                    r2s = []
                    ct+=1
                    for macro in macros:
                        macro = eval(macro)
                        r2 = macroecotools.obs_pred_rsquare(np.log10(SAD),np.log10(macro)) # r-squared for the 1:1 line (obs vs. exp) 
                        r2s.append(r2)
                        
                    density = gaussian_kde(r2s)
                    n = len(r2s)
                    xs = np.linspace(0.0,1.0,n)
                    density.covariance_factor = lambda : .5
                    density._compute_covariance()
                    D = [xs,density(xs)]
                    plt.xlim(0.0, 1.0)
                    plt.setp(ax, xticks=[0.2,0.4,0.6,0.8,1.0])
                    plt.plot(D[0],D[1],lw=0.5,alpha=0.3)
                    data.close()
                else:data.close()
        
        plt.subplots_adjust(wspace=0.3, hspace=0.3)
        print ct,'usable sites in',dataset
        i+=1
    plt.savefig('Figure8.png', dpi=1200, pad_inches=0)

####################################################################################################################
######   Two functions written by Justin Kitzes. One generates uniform random macrostates from the feasible set#####
######   The other function generates a frequency distribution yielding the frequency with which multiplicities ####
######   of an integer appear within the feasible set ##############################################################


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
    
#################################################################################
############# The island of misfit toys #########################################

def next_partition(p):
  if max(p) == 1:
    #return [sum(p)]
    return [0]
  p.sort()
  p.reverse()
  q = [ p[n] for n in range(len(p)) if p[n] > 1 ]
  q[-1] -= 1
  if (p.count(1)+1) % q[-1] == 0:
    return q + [q[-1]]*((p.count(1)+1) / q[-1])
  else:
    return q + [q[-1]]*((p.count(1)+1) / q[-1]) + [(p.count(1)+1) % q[-1]]
