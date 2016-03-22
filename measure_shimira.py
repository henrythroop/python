# Program to calculate square area of Shimira apartment
# Idea is to take some PNG files created in photoshop, and count up the # of white and non-white pixels
# in them. We measure area from those.
#
# HBT 10-Aug-2015

import sys # for sys.sdout.write()

import pdb
import subprocess

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
import sys

import numpy as np

import PIL		 # Python Imaging Library
from PIL import Image

from matplotlib.path import Path
import matplotlib.patches as patches

file_total7  = '/Users/throop/Desktop/Shimira Area 7 Total.png'
file_usable7 = '/Users/throop/Desktop/Shimira Area 7 Usable.png'

file_total8  = '/Users/throop/Desktop/Shimira Area 8 Total.png'
file_usable8 = '/Users/throop/Desktop/Shimira Area 8 Usable.png'

tot7 = np.array(Image.open(file_total7).getdata())	# np.array() is the slow portion of the call here -- don't know why
usa7 = np.array(Image.open(file_usable7).getdata())

tot8 = np.array(Image.open(file_total8).getdata())
usa8 = np.array(Image.open(file_usable8).getdata())

white = [255,255,255,0]

# Now count up the total # of non-white pixels

frac7 = 1. * np.sum(usa7 != white) / np.sum(tot7 != white)
frac8 = 1. * np.sum(usa8 != white) / np.sum(tot8 != white)

a_tot = 51 * 38 # From embassy's drawings. This is the area defined in the 'total' image. I am using this basically
                # to normalize the dimensions.

sf_tot = a_tot * (frac7 + frac8)

print "Floor 7: " + repr(frac7 * a_tot) + ' / ' + repr(a_tot)
print "Floor 8: " + repr(frac8 * a_tot) + ' / ' + repr(a_tot)
print "Total:   " + repr((frac8 + frac7) * a_tot) + ' / ' + repr(a_tot*2)

