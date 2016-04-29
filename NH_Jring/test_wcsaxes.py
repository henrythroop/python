# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 21:01:17 2016

@author: throop
"""

# Short program to test WCSaxes function on a New Horizons LORRI image
# 
# Problem: When plotted, y axis is labeled as -18 50, which is a Dec value. 
#          But, the image is plotted such that the correct y axis should be RA.
# DS9 plots the image properly and displays proper WCS coords.
# Changing wcs.crval offsets the image in the expected directions.
# This suggests that the error is in WCSaxes, rather than the header or astropy.wcs.
#
# Henry Throop 3-Apr-2016

import matplotlib.pyplot as plt
import astropy
import wcsaxes
import numpy as np
from   astropy.io import fits
from   astropy.wcs import WCS

file = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/lor_0034600923_0x630_sci_1.fit'

hdulist = fits.open(file)
header = hdulist[0].header

image = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'

w = WCS(file)                 
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1, projection=w)

# Label the RA and Dec axes.

ra  = ax.coords['ra']
dec = ax.coords['dec']
ra.set_axislabel('RA')
dec.set_axislabel('Dec')

fig1 = plt.imshow(np.log(np.clip(image.data,100,400)), cmap='Greys_r')

fig.savefig('test_wcsaxes.png')

#hdulist.close()