# -*- coding: utf-8 -*-
"Standard HBT function library"
# Created on Wed Nov  5 20:59:35 2014

# To use: 'import HBT as hbt'
# Could do 'from HBT import *', but I don't think importing into global space is a good idea.
#
# Then do 'val = hbt.wheremin(np.arange(100))'

# @author: throop

import numpy as np
import astropy.modeling


# We want to define these as functions, not classes
# They are all general-purpose functions.
# I want to define this all as a class, I guess? Class HBT?
# But when I import it, I want to do this:
#  import HBT as hbt
#  hbt.maxval(arr)    
#########
# Function for wheremin()
#########

def wheremin( arr ):
   "Determines the index at which an array has its minimum value"
   index = np.where(arr == np.amin(arr))
   return (np.array([index]).flatten())[0] # Crazy syntax returns either a scalar, or the 0th element of a vector

def commonOverlapNaive(text1, text2):  
  x = min(len(text1), len(text2))  
  while x > 0:  
    if text1[-x:] == text2[:x]:  
      break  
    x -= 1  
  return x 

def longest_common_substring(S,T):
    "Given two strings, returns the longest common substring as a set"
    "http://www.bogotobogo.com/python/python_longest_common_substring_lcs_algorithm_generalized_suffix_tree.php"
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return lcs_set.pop() # Original function returned lcs_set itself. I don't know why -- I just want to extract one element.

def sfit(arr, degree=3, binning=16): # For efficiency, we downsample the input array before doing the fit.
    "Fit polynomial to a 2D array, aka surface."
    
    shape_small = (np.size(arr,0)/binning, np.size(arr,1)/binning)
    shape_big   = np.shape(arr)

# Create x and y arrays, which we need to pass to the fitting routine

    x_big, y_big = np.mgrid[:shape_big[0], :shape_big[1]]
    x_small = np.resize(x_big, shape_small, order=1, preserve_range=True)
    y_small = np.resize(y_big, shape_small, order=1, preserve_range=True)
    
    arr_small = np.resize(arr, shape_small, order=1, preserve_range=True)
    p_init = astropy.modeling.models.Polynomial2D(degree=degree)

# Define the fitting routine

    fit_p = astropy.modeling.fitting.LevMarLSQFitter()
        
    with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')

# Do the fit itself
        
        poly = fit_p(p_init, x_small, y_small, arr_small)

# Take the returned polynomial, and apply it to our x and y axes to get the final surface fit

    surf_big = poly(x_big, y_big)
    
    return surf_big
    
##########
# Function to determine if an entered string is a number. 
# From http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python
##########

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

