
import sys # for sys.stdout.write
from IPython import embed

from astropy.convolution import convolve, Box1DKernel
from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
import astropy
from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 
import pdb
import cspice
import scipy.ndimage.interpolation

import numpy as np
from   pylab import *
from   scipy.optimize import curve_fit

file_in = '/Users/throop/Data/SALT_Pluto_2015/product/spect_pluto_merged_2015-06-30.txt'
d = loadtxt(file_in, delimiter = ',')

d[d<0] = 0

w = d[:,0] # Wavelength
f = d[:,1] # Flux

num_w = len(w)

num_w2 = 1000
num_w3 = 2000

#   plot(w,f)
#   plt.show()

ff = copy(f)
ff[f == 0] = None

ffc = np.convolve(ff, 1/100. + np.zeros(100), mode = 'same')

plot(w, ffc)
plot(w, ff + 0.1)
plot(w, f + 0.2)
plt.show()







