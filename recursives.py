#!/usr/bin/env sage -python

from sage.all import *
import sys
import os
from os import path, access, R_OK  # W_OK for write permission
import  matplotlib.pyplot as plt
from pylab import *
import numpy as np
import re


""" Here are several recursion functions that accomplish integer partitioning related tasks """

def CountCombinations(N, Minimal): # A recursive function to find the number of partitions with Minimal as the smallest part
  temp = 1
  if N<=1:return 1
  for i in range(1,int(floor(N/2.0))+1):
      if i>=Minimal:
          temp = temp + CountCombinations(N-i, i)
  return temp

def recursive_sum(nested_num_list): # A simple recursive to find the sum of a list of numbers 
    _sum = 0
    for element in nested_num_list:
        if type(element) == type([]):
            _sum = _sum + recursive_sum(element)
        else:
            _sum = _sum + element
    return _sum

def parts(N,largest):   # A recursive function that finds the number of partitions for N with 'largest' as the largest part
    if largest==0: return 0
    if N==0: return 1
    if N<0: return 0
    return parts(N, largest-1) + parts(N-largest,largest)


def parts_NSX(N,S,X):  # A recursive function that finds the number of partitons of N with S parts having X as the largest part
     
     if N==0 and S==0:return 1
     if N<=0 or S<=0 or X<=0:return 0
     if N>0 and S>0 and X>0:
         _sum = 0
         for i in range(0,S+1):
             _sum += parts_NSX(N-i*X, S-i, X-1)
         return _sum    
                  
     
     
     
     
     
     
