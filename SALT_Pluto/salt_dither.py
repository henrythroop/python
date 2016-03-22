# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 12:00:21 2014
# Program reads all of the SALT dataset and makes sense of it.
# Makes plot of which images were taken when
# Creates a post-facto 'observing log' for the SALT data
# Doesn't do any data analysis at all. Does all it work just based on reading & analyzing FITS headers.
#
# HBT Aug-2014

@author: throop
"""

import pdb
import glob
import os.path
from   subprocess import call
import astropy
from   astropy.io import fits
import matplotlib.pyplot as plt # pyplot
import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *  # So I can change plot size: 
import cspice
                     # http://stackoverflow.com/questions/332289/how-do-you-change-the-size-of-figures-drawn-with-matplotlib

year_obs = 2015		# Set this to 2014 or 2015. That will indicate which dataset to analyze

dir_data = "/Users/throop/Dropbox/data/SALT_Pluto_2015/product"

files = [dir_data + '/' + 'mbxgpP201505030051.fits',	# Pluto Blue
         dir_data + '/' + 'mbxgpP201505030053.fits']	# Pluto Red

for file in files:
  hdulist = fits.open(file)

  print 'Number of planes: ' + repr(len(hdulist))


#   jd.append(hdulist[0].header['JD'])
#   hdulist.close()

quit

