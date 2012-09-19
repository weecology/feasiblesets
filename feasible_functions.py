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



########################################################################################################
######   A Section of evenness indices and related functions ###########################################


def simplest_gini(x): #x is a vector of integers
    """This script was obtained from: https://subversion.american.edu/aisaac/notes/blinder_slides.xhtml.
    It yields Gini's coefficient of inequality, a common metric in economics for characterizing inequality
    in distributions of wealth"""
    #initialize yn and ysum
    yn, ysum, countx = 0.0, 0.0, 0
    #compute yN and ysum
    for xn in sorted(x):
      yn = (yn + xn)
      ysum = (ysum + yn)
      countx = (countx + 1)
    #compute area below Lorenz curve
    B = ysum / (countx * yn)
    return(1 - 2*B)

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

def feasible_set_Evar(N,S):
    """ Compute Evar for each macrostate in the feasible set of N and S """
    Evars = []
    for partition in Partitions(N,length=S):
        sad = list(partition)
        Evar = e_var(sad)
        Evars.append(Evar)
    return Evars 

########################################################################################################
######   A Section devoted to finding random macrostates/integer partitions ############################

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
    while ct < 5000:
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
    
def get_rand_sample(NS_combo):
    unique_SADs = []
    pool = Pool()
    unique_SADs = pool.map(worker2, [NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo])
    pool.close()
    pool.join()
    return unique_SADs

def get_macrostates(datasets):
    NS_combos = []
    for dataset in datasets:
        # To prevent the script from tripping on the last line, make the following the last line of the datafile: XXX 111     
        DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
        ct = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        num = 0
        for d in DATA:
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0: SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 9 and sum(SAD) < 100000:
                    NS_combos.append([sum(SAD),len(SAD)])
                    ct+=1
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0: SAD.append(abundance)
                
        print ct,'usable sites in',dataset       
    NS_combos = [list(x) for x in set(tuple(x) for x in NS_combos)]
    NS_combos = sorted(NS_combos,key=itemgetter(0))
    _len = len(NS_combos)
    
    # FINDING MACROSTATES FOR THE ABOVE LIST OF N-S COMBINATIONS
    
    while len(NS_combos) > 0:
        ct = 0
        for NS_combo in NS_combos:    
            ct+=1
            N = int(NS_combo[0])
            S = int(NS_combo[1])
            OUT = open('/home/kenlocey/combined/' + str(N) + '-' + str(S) + '.txt','a+')
            macros = len(set(OUT.readlines()))
            
            if macros < 500 and macros != number_of_partitions(N,S):
                print N,S,ct,_len
                rand_macros = get_rand_sample(NS_combo) # Use multiprocessing
                for i in rand_macros:
                    for p in i:
                        print>>OUT,p
                
                #rand_macros = random_parts(N,S,500)
                #for i in rand_macros:
                #    print>>OUT,i
                OUT.close()
                
            else:OUT.close()
    return

########################################################################################################
##### A Section devoted to examining randoms samples of feasibles sets and empirical data ##############


def get_modal(_list):
    """ Finds the kernel density function across a sample """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    d = 0
    maxd = 0.0
    while d < len(D[1]):
        if D[1][d] > maxd:
            maxd = D[1][d]
            exp_mode = D[0][d] # expected mode is the Simpsons evenness value with the greatest kernel density
        d += 1
    return exp_mode
    
    
def get_hottest_SAD(unique_SADs):
    """ Find the SAD in a random sample with the greatest average commonness 
        among its ranked abundance states. This SAD is taken to represent the 
        central tendency of the set, based on the SAD shape. """
            
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
    a1 = 0
    v1 = 0
    n = 0
    unique_SADs = []
    for p in Partitions(N,length=S):
        unique_SADs.append(list(p))
    for rad in unique_SADs:
        n += 1
        #print n,' ',LU
        C_in_rank = []
        ct1 = 0
        for a in rad: # for each rank
            c = 0
            for x in unique_SADs: 
                if a == x[ct1]:
                    c += 1
            C_in_rank.append(c)
            ct1 += 1
        a2 = mean(C_in_rank)
        log_abs = []
        for a in rad:log_abs.append(ln(RDF(a)))
        v2 = variance(log_abs)  
        if a2 > a1:
            a1 = a2
            v1 = v2
            xRAD = rad
        elif a2 == a1:
            if v2 < v1:
                a1 = a2
                v1 = v2
                xRAD = rad
    return xRAD    

     
def generate_obs_pred_data(datasets):
    
    for dataset in datasets:
        """ To prevent the script from tripping on the last line, make the last line of the datafile any sequence of non-data
            related alphanumeric characters (i.e. 0X0)
        """    
        OUT1 = open('/home/kenlocey/'+dataset+'/'+ dataset + '_obs_pred.txt','w')
        OUT1.close()
        DATA = open('/home/kenlocey/'+dataset+'/'+ dataset + '-data.txt','r')
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
        print len(SADs),'usable sites in',dataset
        
        num = 0
        for site in SADs:
            site_name = site[0]
            N = site[1]
            S = site[2]
            SAD = site[3]
            if max(SAD) > 1:
                
                unique_SADs = []
                PATH = '/home/kenlocey/combined/'+str(N)+'-'+str(S)+'.txt'
                if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                    data = open(PATH,'r')
                    macrostates = data.readlines()
                    for macrostate in macrostates:
                        sad = eval(macrostate)
                        unique_SADs.append(sad)
                    data.close()
                #print N,S,len(unique_SADs)
                if len(unique_SADs) < 100:
                    continue
                if len(unique_SADs) > 500:
                    unique_SADs = random.sample(unique_SADs,500)
                unique_SADs = [list(x) for x in set(tuple(x) for x in unique_SADs)]
                num += 1
                expSAD = get_hottest_SAD(unique_SADs) # The expected SAD from the random sample  
                SAD.sort()
                SAD.reverse()
                
                r2 = macroecotools.obs_pred_rsquare(np.log10(SAD), np.log10(expSAD))
                print dataset,N,S,site_name,r2,' ',berger_parker(SAD),berger_parker(expSAD),' ',simplest_gini(SAD),simplest_gini(expSAD),' ',e_var(SAD),e_var(expSAD) 
                
                ct = 0
                OUT1 = open('/home/kenlocey/'+dataset+'/'+ dataset + '_obs_pred.txt','a')
                while ct < len(expSAD): # write to file, by cite, observed and expected ranked abundances
                    print>>OUT1, site_name,SAD[ct],expSAD[ct]
                    ct += 1
                OUT1.close()
        
        DATA.close()
        print 'FINISHED: ',dataset,' ',num,'\n'  
        

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
    
        

def obs_pred_r2_multi(datasets, data_dir='/home/kenlocey/'): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)' 
    for i, dataset in enumerate(datasets):
        obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt') 
        obs = ((obs_pred_data["obs"]))
        pred = ((obs_pred_data["pred"]))
        print dataset,' ',macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))     


def plot_obs_pred_sad(datasets, data_dir='/home/kenlocey/', radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
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
        
    plt.savefig(dataset + 'obs_pred_plots.png', dpi=400, bbox_inches = 'tight', pad_inches=0)  
    


def get_kdens(N,S):
    """ Finds the kernel density function for Evar across all macrostates 
        of a feasible set based on N and S """
    Evar_obs = feasible_set_Evar(N,S)
    density = gaussian_kde(Evar_obs)
    n = len(Evar_obs)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
    
def get_kdens_obs(SADs):
    """ Finds the kernel density function for Evar across a sample of SADs
        from a feasible set based on N and S """
    Evars = Evars_sample(SADs)
    density = gaussian_kde(Evars)
    n = len(Evars)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

   
def Evar_kdens_full_feasibles(N,S):
    """ Plot kernel density curves of Evar for macrostates of feasible sets
        based on different N and S """
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    D = get_kdens(N,S)
    plt.xlim(0.0, 1.0)
    plt.plot(D[0],D[1],color='black',lw=5)
    D = get_kdens(N,S+10)
    plt.xlim(0.0, 1.0)
    plt.plot(D[0],D[1],color='gray',lw=5)
    plt.axvline(x=0.673,ymin=0,ymax=10,color='black',ls='--',lw=3) # plot a vertical line at the mode
    plt.setp(ax, xticks=[0.2,0.6,1.0],yticks=[0,2,4,6])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    
    plt.savefig('Figure3.png', dpi=400, pad_inches=0)     


def Evar_kdens_feasible_vs_obs(N,S,dataset):
    """ Plot kernel density curves of Evar for macrostates of a feasible set
        and a sample of macrostates """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
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
                SADs.append(SAD)
                
            SAD = []
            abundance = int(re.findall(r'\d*\S$',d)[0])
            if abundance > 0:SAD.append(abundance)
    print len(SADs),'usable sites in',dataset
    DATA.close()
    D = get_kdens(N,S)
    plt.plot(D[0],D[1],color='black',lw=5)
    D = get_kdens_obs(SADs)
    plt.plot(D[0],D[1],color='gray',lw=5)
    
    plt.setp(ax, xticks=[0.2,0.6,1.0],yticks=[1,3,5])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    plt.savefig('Figure2.png', dpi=400, pad_inches=0) 


def get_NS_combos(datasets):
    
    NS_combos = []
    for dataset in datasets:
        
        DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
        ct1 = 0
        ct2 = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        num = 0
        for d in DATA:
            ct1+=1
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 9:
                    NS_combos.append([sum(SAD),len(SAD)])
                    ct2+=1
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
                
    return NS_combos
    

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


def get_all_SADs(N,S): # Figure 1 Locey and White (20??) 

    fig = plt.figure()
    ax = fig.add_subplot(1,3,1)
    for p in Partitions(N,length=S):
        sad_data = list(p)
        (bins, n) = histOutline(sad_data)
        plt.plot(bins, n, 'k-',lw=1,alpha=0.02)
    
    ct = 0
    while ct == 0:
        sad_data = Partitions(N).random_element()
        if len(sad_data) == S:
            ct+=1
            (bins, n) = histOutline(sad_data)
            plt.plot(bins, n, 'r-',lw=2)
    plt.xlim(1,N-S)
    plt.ylim(0,S)
    plt.setp(ax, xticks=[10,20,30,40],yticks=[0,2,4,6,8,10])
    print '1st subplot finished'
    ax = fig.add_subplot(1,3,2)
    for p in Partitions(N,length=S):
        RAD = list(p)
        plt.plot(RAD,color='black',lw=1,alpha=0.02)
    ct = 0
    while ct == 0:
        p = Partitions(N).random_element()
        if len(p) == S:
            ct+=1
            RAD = list(p)
            plt.plot(RAD,color='red',lw=2)    
    plt.ylim(0,40)
    plt.setp(ax, xticks=[0,4,9])
    locs,labels = xticks()
    xticks(locs, map(lambda x: int(x+1), locs))
    plt.setp(ax,yticks=[10,20,30,40])
    print '2nd subplot finished'
    ax = fig.add_subplot(1,3,3)
    D = get_kdens(N,S)
    plt.plot(D[0],D[1],color='black',lw=3)
    plt.setp(ax, xticks=[0.2,0.6,1.0],yticks=[1,2,3])
    
    print 'done'
    plt.subplots_adjust(bottom=0.3, top=0.6,wspace=0.5)
    plt.savefig('Fig1-'+str(N)+'-'+str(S)+'.png', dpi=400, pad_inches=0)
    
    
def get_SADs(dataset):

    DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
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
            if len(SAD) > 9:
                SAD.sort()
                SAD.reverse()
                SADs.append(SAD)
                #if len(SADs) >= 400:break
                ct2+=1
            SAD = []
            abundance = int(re.findall(r'\d*\S$',d)[0])
            if abundance > 0:SAD.append(abundance)
    DATA.close()
    return(SADs)
    
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
    i = 0
    for dataset in datasets:
        ax = fig.add_subplot(3,3,i+1)
        print dataset
        SADs = []
        
        x_list_evar = []
        x_list_NHC = []
        x_list_EQ = []
        y_list_evar = []
        y_list_NHC = []
        y_list_EQ = []
        
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
                        macros = random.sample(macros,500)
                    obs_evar = e_var(SAD)
                    x_list_evar.append(obs_evar)
                    obs_EQ = EQ_evenness(SAD)
                    x_list_EQ.append(obs_EQ)
                    obs_NHC = NHC_evenness(SAD)
                    x_list_NHC.append(obs_NHC)
                    macrostates = []
                    for macro in macros:
                        macrostates.append(eval(macro))
                    sample_evar = []
                    sample_EQ = []
                    sample_NHC = []
                    for macro in macrostates:
                        sample_EQ.append(EQ_evenness(macro))
                        sample_evar.append(e_var(macro))
                        sample_NHC.append(NHC_evenness(macro))
                    
                    percentile_evar = stats.percentileofscore(sample_evar,obs_evar)
                    percentile_E = stats.percentileofscore(sample_EQ,obs_EQ)
                    percentile_NHC = stats.percentileofscore(sample_NHC,obs_NHC)
                    y_list_evar.append(percentile_evar)
                    y_list_EQ.append(percentile_E)
                    y_list_NHC.append(percentile_NHC)
                    data.close()
                    #if len(x_list_evar) >= 20:break
                else:data.close()
        
        slope,intercept,r_value,p_value,std_err = stats.linregress(x_list_EQ,y_list_EQ)
        print 'percentile EQ: r-value:',r_value,'p-value:',p_value,'slope:',slope,' r-squared:',r_value**2
        #m,b = np.polyfit(x_list_EQ,y_list_EQ,1)
        #plt.plot(x_list_EQ, np.array(x_list_EQ)*m +b, c='b',lw=3) 
        plt.scatter(x_list_EQ,y_list_EQ,c='b',marker='o',s=10,lw=0,alpha=0.3)
        
        slope,intercept,r_value,p_value,std_err = stats.linregress(x_list_evar,y_list_evar)
        print 'obs-pred Evar: r-value:',r_value,'p-value:',p_value,'slope:',slope,' r-squared:',r_value**2
        #m,b = np.polyfit(x_list_evar,y_list_evar,1)
        #plt.plot(x_list_evar, np.array(x_list_evar)*m +b, c='r',lw=3) 
        plt.scatter(x_list_evar,y_list_evar,c='r',marker='o',s=10,lw=0,alpha=0.3)       
        
        plt.xlim(0,1)
        plt.ylim(0,100)
        plt.setp(ax, xticks = [0.2,0.4,0.6,0.8], yticks=[20,40,60,80])
        plt.tick_params(axis='both', which='major', labelsize=8)
        
        # Create inset for kdens of site level comparisons of evenness to the feasible set
        slope,intercept,r_value,p_value,std_err = stats.linregress(x_list_NHC,y_list_NHC)
        print 'percentile NHC: r-value:',r_value,'p-value:',p_value,'slope:',slope,' r-squared:',r_value**2,'\n'
        axins = inset_axes(ax, width="40%", height="40%", loc=2)
        plt.scatter(x_list_NHC, y_list_NHC, c='0.25',marker='o',lw=0,s=10,alpha=0.3)
        plt.setp(axins, xticks=[], yticks=[])
        
        i+=1
        
    plt.savefig('percentile.png', dpi=400, bbox_inches = 'tight', pad_inches=0)  
    
    
    
def get_500_RADs(NS_combos): # FIGURE 4 LOCEY AND WHITE (20??)
    
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



def get_common_combos(datasets):
    n = 0
    all_NS_combos = get_NS_combos(datasets)
    unique_NS_combos = [list(x) for x in set(tuple(x) for x in all_NS_combos)]
    common_combos = []
    for combo in unique_NS_combos:
        if all_NS_combos.count(combo) >= 2 and combo[0] >= 100:
            num = all_NS_combos.count(combo)
            print combo,num
            n+=num            
            common_combos.append(combo)
    print 'number of unique N-S combos:',len(common_combos),'  number of sites:',n
    
    
def FIA_NS_combos(NS_combos,dataset):
    """ Plot kernel density curves of species evenness for macrostates of feasible sets
        based on different N and S """
    
    i = 1
    fig = plt.figure()
    for combo in NS_combos:
        ax = fig.add_subplot(4,4,i)
        N = int(combo[0])
        S = int(combo[1])
        DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
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
        print len(SADs),'usable sites for N =',N,'and S =',S
        DATA.close()
        D = get_kdens(N,S)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='black',lw=3)
        maxd = 0.0
        xspot = 0
        d = 0
        while d < len(D[1]):
            if D[1][d] > maxd:
                maxd = D[1][d]
                xspot = D[0][d]
            d += 1
        plt.axvline(x=xspot,ymin=0,ymax=10,color='black',ls='--',lw=1.5) # plot a vertical line at the mode
        
        D = get_kdens_obs(SADs)
        plt.xlim(0.0, 1.0)
        plt.plot(D[0],D[1],color='gray',lw=3)
        plt.setp(ax, xticks=[0.2,0.8])
        plt.subplots_adjust(wspace=0.35, hspace=0.35)
        i+=1  
    plt.savefig('FIA_feasible_vs_obs.png', dpi=400, pad_inches=0)     
    

def pairwise_r2_obs_feasible(datasets):
    
    i = 1
    fig = plt.figure()
    for dataset in datasets:
        ax = fig.add_subplot(3,3,i)
        P_errs = []
        DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
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
                if len(SAD) > 9:
                    SAD.sort()
                    SAD.reverse()
                    SADs.append(SAD)
                #if len(SADs) > 50:break
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
        DATA.close()
        ct = 0

        for SAD in SADs:
            N = sum(SAD)
            S = len(SAD)
            PATH = '/home/kenlocey/combined/' + str(N) + '-' + str(S) + '.txt'
            if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
                data = open(PATH,'r')
                macros = set(data.readlines())
                if len(macros) >= 100:
                    if len(macros) >= 500:
                        macros = random.sample(macros,500)
                    r2s = []
                    perc_errs = []
                    ct+=1
                    for macro in macros:
                        macro = eval(macro)
                        r2 = macroecotools.obs_pred_rsquare(np.log10(SAD),np.log10(macro)) # r-squared for the 1:1 line (obs vs. exp) 
                        r2s.append(r2)
                        #obsEvar = e_var(SAD)
                        #expEvar = e_var(macro)
                        #percent_err = 100*(abs(obsEvar - expEvar)/obsEvar) 
                        #perc_errs.append(percent_err)
                        #if len(perc_errs) >= 10:break
                    #P_errs.append(perc_errs)    
                    density = gaussian_kde(r2s)
                    n = len(r2s)
                    xs = np.linspace(0.0,1.0,n)
                    density.covariance_factor = lambda : .5
                    density._compute_covariance()
                    D = [xs,density(xs)]
                    plt.xlim(0.0, 1.0)
                    plt.setp(ax, xticks=[0.2,0.4,0.6,0.8,1.0])
                    plt.plot(D[0],D[1],color='black',lw=0.5,alpha=0.3)
                    data.close()
                else:data.close()
        
        plt.subplots_adjust(wspace=0.3, hspace=0.3)
        print ct,'usable sites in',dataset
        i+=1
    plt.savefig('Figure8.png', dpi=400, pad_inches=0)
    
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
    

        
def compare(size,datasets):
    
    NS_combos = []
    for dataset in datasets:
        
        DATA = open('/home/kenlocey/' + dataset + '/' + dataset + '-data.txt','r')
        ct1 = 0
        ct2 = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group()
        
        m2 = int(re.findall(r'\d*\S$',d)[0])
        SAD = [int(m2)]
        num = 0
        for d in DATA: # for each line in the dataset
            ct1+=1
            m1 = re.match(r'\A\S*',d).group()
        
            if m1 == m0:
                m2 = int(re.findall(r'\d*\S$',d)[0])
                if m2 > 0:SAD.append(m2)
            
            else:
                site_name = m0
                m0 = m1
                if len(SAD) > 9:
                    NS_combos.append([sum(SAD),len(SAD)])
                    ct2+=1
                SAD = []
                abundance = int(re.findall(r'\d*\S$',d)[0])
                if abundance > 0:SAD.append(abundance)
    
        NS_combos = [list(x) for x in set(tuple(x) for x in NS_combos)]
    
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
            
            D = get_kdens_obs(SADs)
            plt.xlim(0.0, 1.0)
            plt.plot(D[0],D[1],color='black',lw=2)    
            
            SADs = []
            SADs = random_parts(N,S,size)
            D = get_kdens_obs(SADs)
            plt.xlim(0.0, 1.0)
            plt.plot(D[0],D[1],color='red',lw=2)
            data.close()
            
            print 'finished:'+str(N)+' '+str(S)
            plt.savefig('/home/kenlocey/ZPICS/N'+str(N)+'-S'+str(S)+'-'+str(size)+'macros.png', dpi=400, pad_inches=0)
        else: print 'no file for '+str(combo[0])+'-'+str(combo[1])

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