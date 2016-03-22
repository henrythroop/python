# Simple program to fix the headers on one set of SALT images.
# As per Steve Crawford e-mail 17-Aug-2015
#
# HBT Sep-2015

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

dir_data = "/Users/throop/Dropbox/data/SALT_Pluto_2015/product"

# Set the file

files = np.array(['mbxgpP201506290032.fits'])		# 31, 32

# Read the file

for file in files:
    hdulist = fits.open(dir_data + '/' + file)

    print 'MASKTYP = ' + hdulist[0].header['MASKTYP']
    print 'MASKID = '  + hdulist[0].header['MASKID']

# Change the header fields. These two are incorrectly set in this file alone. All other files are OK, and
# all other fields are OK.

    hdulist[0].header['MASKTYP'] = 'LONGSLIT          '
    hdulist[0].header['MASKID'] = 'PL0200N001        '

    file_out = dir_data + '/' + '0032.fits'

# Write file back to disk (under new name)

    hdulist.writeto(file_out)
    print 'Wrote: ' + file_out

    hdulist.close()
