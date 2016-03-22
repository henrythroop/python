# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:10:23 2014

@author: throop
"""

import pdb
import glob
import os.path
from   subprocess import call
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt # pyplot
import numpy as np
from pandas import DataFrame, Series
import pandas as pd
from pylab import *  # So I can change plot size: 
import png

file_in = '/Users/throop/Dropbox/CCLDAS/Run2 Oct13/SEM Roy 2014/1-E5_Pan.tif'
file_in = '/Users/throop/Desktop/1-E5_Pan.tif'

png.reader(file_in)

from libtiff import tiff # pip install libtiff

# to open a tiff file for reading:
file_in = '/Users/throop/CCLDAS/Run2 Oct13/SEM Roy 2014/Throop 1-E5 23Sep14/1-E5_Pan.tif'

tif = tiff.open(file_in, mode='r') # Use libtiff, which is a python wrapper to the C library I installed.
tif = TIFFfile(file_in) # Use native PyTIFF


# to read an image in the currect TIFF directory and return it as numpy array:

image = tiff.read_image()


p=png.Reader(file_in)
pp=p.read(lenient=True)
# FormatError: FormatError: Chunk �C� is too large: 3821738563.
