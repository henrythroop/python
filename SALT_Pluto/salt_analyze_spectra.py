# This routine will read in the fully processed and merged / combined SALT spectra.
# It will then do science on them: measure bands, do smoothing, etc.
#
#
# HBT 29-Oct-2014
# HBT  2-Nov-2015 Updated for 2015 data
# HBT    Dec-2015 Rewrote significantly to work with the 2014 and 2015 data together. It should be used, not the old one.

import pdb
import sys # for sys.stdout.write

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
# import astropy
# from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 
from   subprocess import call   # Shell output
from   subprocess import check_output   # Shell output

import numpy as np
# from   pandas import DataFrame, Series
from   pylab import *
#from   astropy.convolution import convolve, Box1DKernel
import cspice

from scipy.optimize import curve_fit

# Wrapper to plot() to extract specified elements. *args and **kwargs are same as IDL __extra .

def plot_elements(data_x, data_y, elements, plot_all = False, *args, **kwargs):
    """ Plot individual elements of an array """
    if (plot_all):                          # plot_all: flag to ignore setting of 'elements', and plot everything
        return plot(data_x, data_y, *args, **kwargs)
    else:
        return plot(data_x[elements], data_y[elements], *args, **kwargs)
    
##########
# Simplest method of fitting a band: 
# Draw a line between two wavelengths surrounding the CH4 band.
# Then get the area under that curve.
##########

def calc_area_curve(data_x, data_y, alam, binning):
    """ Calculate area between a straight line connecting two points, and a spectrum """
    
    bin_0  = np.min(np.where(data_x > alam[0]))
    bin_1  = np.max(np.where(data_x < alam[1])) 
    
    data_y_smooth = smooth_boxcar(data_y, binning)
    
    data_x_use = data_x[bin_0:bin_1]
    data_y_use = data_y_smooth[bin_0:bin_1]

    dy = data_y_use[-1] - data_y_use[0]		# Get the overall dx and dy from one end of the band to the other
    dx = data_x_use[-1] - data_x_use[0]

    m = dy / dx		          # Slope
    b = data_y_use[0] - m * data_x_use[0] # Y-intercept

    line = m * data_x_use + b			# Get an equation for the line, so we can plot it.
    						# We don't fit to the line or anything.

    area = sum(line - data_y_use) 		# Sum the total area. Line - data because data is below. Don't use ABS().
    
    line_out = m * data_x + b    # A fit to the line, covering the whole input data range
    
    print 'calc_area_curve: x0 = ' + repr(alam[0]) + ', x1 = ' + repr(alam[1]) + ', y0 = ' + repr(line_out[bin_0])
    return (area, line_out)   
    
#################    
# Returns the 5 parameter (or 3-parameter) best fit
# NB: data_x should usually be the flux *ratio*

def calc_fit_gaus(data_x, data_y, alam, initial, constrained=True, binning=1 ):
    
    (m, a, x0, b, sigma) = initial # Import the initial guesses. Note that we re-calculate m and b here in this function.
    
    is_good = (data_y > 0)	# For this observation, ignore any wavelengths where flux is zero

    data_x_full = data_x[is_good]		# data_{xy}_good are where there is non-zero flux

    data_y_full = smooth_boxcar(data_y[is_good], binning)

# Remove the Inf's from the dataset. This is easier working with them than it is in IDL!

    data_x_good = data_x_full[data_y_full != Inf]	# data_{xy}_real is where there is non-zero *and* non-NaN flux
    data_y_good = data_y_full[data_y_full != Inf]

    bin_0     = np.min(np.where(data_x_good > alam[0]))
    bin_1     = np.max(np.where(data_x_good < alam[1]))

    data_y_short = data_y_good[bin_0:bin_1]
    data_x_short = data_x_good[bin_0:bin_1]
    
# First get a good guess for the slope
# Do this by crudely taking the values and start and end of dataset.
# No need for accuracy since it'll be fit later 

    dy = data_y_short[-1] - data_y_short[0]
    dx = data_x_short[-1] - data_x_short[0]

    m = dy / dx		      # Slope
    b = data_y_short[0] - m * data_x_short[0] # Y-intercept

    line = m * data_x_short + b

# Stuff the guesses into an array

    p0_7200             = [m, a, x0, b, sigma]
    p0_7200_constrained = [m,        a, b    ]
  
# Now do the fits.

# If we want constrained (that is, three parameters)
 
    if (constrained == True): 
        try:
            popt,pcov = curve_fit(gaus_constrained, data_x_short, data_y_short, p0=p0_7200_constrained)  # Output is (optimized, covariance)
            success = True
            fits = popt
            (m, a, b) = popt
            data_fit_out = gaus_constrained(data_x, m, a, b)
        except RuntimeError:
            success = False
            data_fit_out = nan + data_x
            fits = EMPTY

# If we want unconstrained (that is, five parameters)

    if (constrained == False):  # 
        try:
            popt,pcov = curve_fit(gaus, data_x_short, data_y_short, p0=p0_7200)  # Output is (optimized, covariance)
            success = True
            fits = popt
            (m, a, x0, b, sigma) = popt
            data_fit_out = gaus(data_x, m, a, x0, b, sigma)
        except RuntimeError:
            success = False
            data_fit_out = nan + data_x
            fits = EMPTY   
    
    return (fits, success, data_fit_out)

#########
# Return common elements of an array. Order of list1 is maintained -- only that elements that are not in list2 are excised.
# Not the fastest way to do it but works well for me since I need to maintain order
# http://stackoverflow.com/questions/2864842/common-elements-comparison-between-2-lists
#########
def common_elements(list1, list2):
    return [element for element in list1 if element in list2]

#########
# Function for converting array to string: [100,200] -> "100 .. 200"
#########

def range2str(arr):
  return repr(arr[0]) + ' .. ' + repr(arr[-1])

##########
# Function to do boxcar smoothing
# Array returned has same length as input array.
# In theory this works basically the same as astropy.convolution.kernels.convolve with Box1DKernel, 
# but that kept giving me python error, so better to bring it in house.
##########

def smooth_boxcar(ydata, binning):
   smoothed = np.convolve(ydata, 1./binning + np.zeros(binning), mode = 'same')
   return smoothed
             
#########
# Function for wheremin()
#########

def wheremin( arr ):
    "Determines the index at which an array has its minimum value"
    index = where(arr == amin(arr))
    return index[0]

##########
# Gaussian function. Used for curve-fitting of spectral lines
##########

# b      Vertical offset at x=0 (that is, y-intercept, aka y0)
# m      Linear slope
# a      Height of Gaussian (can be negative)
# x0     Position of Gaussian
# sigma  Width of Gaussian

def gaus(x,m,a,x0,b,sigma):
  return b + m*x + a*exp(-(x-x0)**2/(2*sigma**2))

def gaus_constrained(x,m,a,b):				# As per Grundy's suggestions, fit a 'constrained' gaussian, which has fewer free parameters.
  x0    = 7.28636414e+03				# Fix the center
  sigma = 2.65981811e+01				# Fix the width
  return b + m*x + a*exp(-(x-x0)**2/(2*sigma**2))

##########
# Function to align_wavelengths()
##########
   
def align_wavelengths( wavelengths, fluxes, wavelength_range, wavelength_center, smoothing=1):
    "Take a bunch of spectra, and slide them right/left until they align"
    
# All spectra share the same identical wavelength axis
# Procedure: 
#   Loop over all wavelengths
#   Extract a portion between min and max
#   If requested, smooth it by 'smoothing' with boxcar    
#   In that spectrum, calc the minimum value between min and max set in wavelength_range
#   Shift it so that position of minimum is set to the wavelength desired
#   ** During the shift, no data are lost. They roll from one end to the other.
#    
# NB: Python parameters are all passed by reference. So the original array gets changed, as intended.   

    bin_start = wheremin(abs(wavelengths - wavelength_range[0]))[0] # Index for starting wavelength in region 
    bin_end   = wheremin(abs(wavelengths - wavelength_range[1]))[0] # Index for ending wavelenngths
    wavel_extract = wavelengths[bin_start : bin_end]      
    index_right = wheremin(abs(wavel_extract - wavelength_center))[0] # Get the proper position

    for i in range(len(fluxes)):   
        fluxes_extract = fluxes[i][bin_start : bin_end].copy()               # Extract the small region to examine
        fluxes_extract_smooth = smooth_boxcar(fluxes_extract, smoothing)
#        fluxes_extract_smooth = roll(fluxes_extract_smooth, -(smoothing-1)/2) # Weird. Sometimes the output array is dfft size from input.
        hw = (smoothing-1)/2 + 1 # halfwidth of the smoothing
        index_wrong = wheremin(fluxes_extract_smooth[hw:-hw])[0] + hw # Get the central position 
        droll = index_right - index_wrong
#        print "index_right = {0}, index_wrong = {1}".format(index_right, index_wrong)
#        print "droll = " + repr(droll)
#        quit       
#        print 'Rolling spectrum #' + repr(i) + ' by ' + repr(droll) + " to match " + repr(wavelength_center)
        fluxes[i] = roll(fluxes[i], droll)
   
    return # No return arguments -- everything is by reference
        
##########
# Start of main program
##########
   
# Set the directory for the data and all analysis. 

EMPTY	= nan   # I did have this set to 999, but I think nan is better. Use it for what it's for. 999 gets parsed as a number.

# The 'set' file, listing combination of Pluto and HD to use.
# Note that at this time, the 'set' file is of fixed length... must be 22 entries, even if we actually 
# don't want to plot all of those.

filename_sets = '/Users/throop/python/salt_interact_settings/salt_2014_2015_11Jan16'

dir_data_2015 = "/Users/throop/Data/SALT_Pluto_2015/product/" 

dir_data_2014 = "/Users/throop/Data/SALT_Pluto_2014/product/" 

#files_hd = glob.glob(dir_data + '/spect_hd*txt')
#files_pluto = glob.glob(dir_data + '/spect_pluto*txt')

# Initialize variables

flux_hd       = []
flux_pluto    = []
wavelength    = [] # I think that all of the files should have an identical wavelength now.
date_pluto    = []	# Date of the Pluto obs. Note that this is always the date of the *evening*... even if it is past midnight
date_hd       = []	# Date of HD. Same evening convention as for Pluto
files_hd   = []
files_pluto= []
year_hd       = []
year_pluto    = []

f = open(filename_sets)
lines = f.readlines()
f.close()

header_comment = (lines[0])[2:-1] # remove the hashtag and newline
header_columns = (lines[1])[:-1]

# Parse the header row

(index, sa, sb, num_lines, num_sets_listboxes, num_sets_optionlists) = header_columns.split(',')

num_lines = int(num_lines)

lines = lines[2:]

# In the particular case of one SALT file I use a lot, the last two lines are duplicates. SALT_INTERACT2.PY has no way to 
# exclude these, so kill them here.

if (num_lines == 22) and (filename_sets == '/Users/throop/python/salt_interact_settings/salt_2014_2015_11Jan16'):
  lines = lines[0:20]
  num_lines = len(lines)

for i in range(num_lines):
#    print "Unpacking line: " + lines[i]
    lines[i].strip()    # Remove whitespace, \n, etc.
    (num, date_pl_i, date_hd_i) = lines[i].strip().split(',')
    date_hd_i = date_hd_i.strip()
    date_pluto_i = date_pl_i.strip()

    date_pluto.append(date_pluto_i)
    date_hd.append(date_hd_i)
    year_pluto_i = int(date_pluto_i[0:4])
    year_hd_i    = int(date_hd_i[0:4])
    year_pluto.append(year_pluto_i)
    year_hd.append(year_hd_i)
    
    files_hd.append('/Users/throop/Data/SALT_Pluto_' + repr(year_pluto[i]) + '/product/spect_hd_merged_' + date_hd_i + '.txt')
    files_pluto.append('/Users/throop/Data/SALT_Pluto_' + repr(year_pluto[i]) + '/product/spect_pluto_merged_' + date_pluto_i + '.txt')
                
year_obs = ''   # In previous version, this was 2014 or 2015. I am keeping it here for compatability, but a blank string

# Load data from all the HD files

for file in files_hd: 
    print "Loading file " + file
    d = loadtxt(file, delimiter=',')
    flux_hd.append(d[:,1])
    l = len(file)
#    date_hd.append(file[l-14:l-4])		# Extract the date string
    wavelength.append(d[:,0])

wavelength = wavelength[0]

# Load data from all the Pluto files
    
for file in files_pluto: 
    print "Loading file " + file
    d = loadtxt(file, delimiter=',')
    flux_pluto.append(d[:,1])
    l = len(file)
#    date_pluto.append(file[l-14:l-4])		# Extract the date string

# Convert to numpy arrays

date_pluto = np.array(date_pluto)
date_hd    = np.array(date_hd)

# Set up plotting parameters

rcParams['figure.figsize'] = 20, 10

colors = ['red', 'blue', 'green', 'purple', 'orange', 'black', 'darkblue', 'salmon', 'pink', 'brown', 'olive', 'violet', 'brown', 'tan', 
          'darkolivegreen', 'pink', 'cyan', 'lightgrey', 'mediumpurple', 'tan', 'sienna', 
                       'rosybrown', 'darkred', 'coral', 'thistle', 'navy', 'yellowgreen', 'darkkhaki',
                       'tomato', 'salmon', 'gold', 'lightblue', 'seagreen']

# Make a plot of all HD obs, so we can pull out the best

xlim = (3500, 9500)
ylim = ((0, 1.2))

#if (year_obs == 2014):
#  index_good_spectra_hd = {0, 1, 2, 3, 4} # 3 looks the best (26-Aug)
#if (year_obs == 2015):
#  index_good_spectra_hd = {0, 1, 2, 3, 4, 5, 6, 7, 8} 
  
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1, color=colors[mod(i,len(colors))])
      plt.text(amax(wavelength), 
               np.median(flux_hd[i][-100:-1]) + i*0.1, 
               '  ' + date_hd[i])

plt.title('HD, SALT ' + repr(year_obs), fontsize=24)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

show()

##########
# Flag some nights as bad. We can then easily exclude these from the plots
##########

is_pluto_good = (date_pluto != '2020')
w_is_pluto_good = where(is_pluto_good)[0]

##########
# Make a plot of region around 0.73 um band
##########

# ** This is the important one. Check RSS Simulator to see where it puts the gap.
# PG0900, actual tilt in FITS files = 20.375
# A: RSS simulator says gap should be at 7150 .. 7200 (for tilt 20.375)
# Data put the gap at                    7175 .. 7225 (for tilt 20.375)
# So it's off by about a half a gap width. That is quite unfortunate. But it's more my fault than SALT's.

# If I would have requested a tilt angle of 40 deg (20 deg), then the gap would be clearly off of the CH4 band, at 7000 .. 7050.
# That is clearly what I intended from my proposal. Fuck.
 
xlim = (6700, 7500)

binning = 10

for i in range(len(files_pluto)):
      plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i] + i*0.1, binning), color=colors[mod(i, len(colors))])
      plt.text(amax(wavelength), 
               np.median(flux_pluto[i][-100:-1]) + i*0.1, 
               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT. Band = CH$_4$, 0.71 $\mu$m - 0.74 $\mu$m', fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim((0.7,2.5))
show()

##########
# Make plot around solar Ha (H alpha) line at 6563 A
##########

xlim = (6520, 6620)
ylim = (0.4, 3.0)
for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[mod(i, len(colors))])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[mod(i, len(colors))])      
      
plt.title("Pluto (bottom) and HD (top). Band H$\\alpha$ = 6563 $\AA$", fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make plot around Cruikshank's 0.62 um CH4 band
##########

if (False):
    xlim = (6100, 6300)
    ylim = (0.4, 3.0)
    for i in range(len(files_pluto)):
	  plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
	  
    for i in range(len(files_hd)):
	  plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
	  
    plt.title("Pluto (bottom) and HD (top). Band = 6200 $\AA$ CH$_4$", fontsize=24)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    show()

##########
# Now center on Ha, and plot again
##########

#wavelengths = wavelength
#fluxes = flux_pluto
#wavelength_range = (6540, 6580)
#wavelength_center = 6563
#smoothing = 10
#
#stop

if (False):
    align_wavelengths(wavelength, flux_pluto, (6540, 6580), 6563, smoothing=10)
    align_wavelengths(wavelength, flux_hd,    (6540, 6580), 6563, smoothing=10)

    xlim = (6520, 6620)
    ylim = (0.4, 3.0)
    for i in range(len(files_pluto)):
	  plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
	  
    for i in range(len(files_hd)):
	  plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
    plt.title("Pluto (bottom) and HD (top). Band H$\\alpha$ = 6563 $\AA$, aligned @ 6563", fontsize=24)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    show()

#stop

##########
# Make plot around Telluric O2 line at 7620 A. See Grundy 1996 @ 332. Used for wavelength calibration.
##########

if (False):

    align_wavelengths(wavelength, flux_pluto, (7580, 7620), 7600, smoothing=10)
    align_wavelengths(wavelength, flux_hd,    (7580, 7620), 7600, smoothing=10)

    xlim = (7570, 7670)
    ylim = (0.2, 2.6)
    binning = 1

    for i in range(len(files_pluto)):
	  plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
	  
    for i in range(len(files_hd)):
	  plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
	  
    plt.title('Pluto (bottom) and HD (top). Band O$_2$ = 7620 $\AA$, aligned @ 7600', fontsize=24)
    plt.xlim(xlim)
    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    plt.ylim(ylim)
    show()
       
#stop
      
      
##########
# Make plot around 0.89 um. HD and Pluto separately. Used for wavelength calibration.
##########

if (False):
    align_wavelengths(wavelength, flux_pluto, (8530, 8550), 8540, smoothing=10)
    align_wavelengths(wavelength, flux_hd,    (8530, 8550), 8540, smoothing=10)

    xlim = (8300, 9300)

    ylim = (0.2, 2.6)
    binning = 1

    for i in range(len(files_pluto)):
	  plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
	  
    for i in range(len(files_hd)):
	  plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
	  
    plt.title('Pluto (bottom) and HD (top). Band O$_2$ = 7620 $\AA$, aligned @ 8540', fontsize=24)
    plt.xlim(xlim)
    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    plt.ylim(ylim)
    show()


##########
# Make a plot of region around 0.89 um band. Ratio.
##########
              
xlim = (8300, 9300)
ylim = (1.0, 4)
binning = 5

index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i] + i*0.2, binning), color=colors[i])
#      plt.text(amax(wavelength), 
#               np.median(flux_pluto[i][-200:-50]) + i*0.2, 
#               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT, 8500 - 8700 and 8700 - 8900', fontsize=24)
plt.xlim(xlim)
plt.ylim((ylim))

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make a plot of region around 0.4 um : Pluto
##########
              
xlim = (3800, 4500)
ylim = (0.0, 0.9)
binning = 5

index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, smooth_boxcar(flux_pluto[i] + i*0.05, binning), color=colors[i])

plt.title('Pluto, SALT, 3500-4000', fontsize=24)
plt.xlim(xlim)
plt.ylim((ylim))

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make a plot of region around 0.4 um : HD
##########

if (False):
    xlim = (3800, 4500)
    ylim = (0.0, 0.9)
    binning = 5

    index_good_spectra_pluto = {0, 1, 2, 3, 4} # 5 is messed up; 11 is noisy
    for i in range(len(files_pluto)):
      if i in index_good_spectra_pluto:
	  plot(wavelength, smooth_boxcar(flux_hd[i] + i*0.05, binning), color=colors[i])
    #      plt.text(amax(wavelength), 
    #               np.median(flux_pluto[i][-200:-50]) + i*0.2, 
    #               '  ' + date_pluto[i])

    plt.title('HD, SALT, 3500-4000', fontsize=24)
    plt.xlim(xlim)
    plt.ylim((ylim))

    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    show()

##########
# Make a plot of region around 0.89 um band, BINNED
##########

if (False):
    align_wavelengths(wavelength, flux_pluto, (8530, 8550), 8540, smoothing=10)
    align_wavelengths(wavelength, flux_hd,    (8530, 8550), 8540, smoothing=10)        

    xlim = (8300, 9300)
    ylim = (0.9, 4)
    binning = 10

    for i in range(len(files_pluto)):
      if i in index_good_spectra_pluto:
          plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[3] + i*0.2, binning), color=colors[i])
          plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])

    plt.title('Pluto/HD, SALT, 8500-8700 and 8700-8900, binned x ' + repr(binning), fontsize=24)
    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)

    file_out = 'spect_ratio_8900_bin{0}.png'.format(binning)
    plt.savefig(file_out)
    print "Wrote: " + file_out
  
    show()

##########
# Look up sub-solar longitudes of Pluto from file
##########

# Read in the conversion table between JD and subsol longitude.
# This file has a complex format (human-readable), so rather than read it in, we just grep it properly
# The subsollon from the file is the mean of all the Pluto obs for that date.
# It is thus more accurate than the one I compute from SPICE in the current file. They are within a few degrees
# of each other.
# 
# Sample line ('#' is not in original)
#
#    JD: 2456830.0, 2014-06-21, night 2, longitude_mean_pl_rss = 333.63548826680176

subsollon     = np.zeros(len(files_pluto))
subsollon_rhr = np.zeros(len(files_pluto))
jd            = np.zeros(len(files_pluto))

for i in range(len(files_pluto)):
    
  file_subsollon = '/Users/throop/Data/SALT_Pluto_' + repr(year_pluto[i]) + '/product/' + 'SALT_Pluto_Observing_Log.txt'

  file = files_pluto[i]

  if (year_obs == '2014'):
    file = file.replace("07-02", "07-01")			# Ugh: the spectra are written with the time after midnight, and the 
    file = file.replace("07-03", "07-02")			# observing log has the time before midnight... I think

  datestr = file[string.find(file, 'merged_') + 7:-4]   # extract just the '2014-08-22' portion of the filename

  cmd = ["grep", datestr, file_subsollon]		# Grep through the Observing Log to find that date
  out = check_output(cmd)
  words = out.split(" ")
  subsollon[i] = float(words[-1])			# Extract last word from line that matches.
                                               # That is the mean Pluto Sub-solar longitude for that date
  subsollon[i] = round(subsollon[i],2)        # Truncate it for clarity
  
  if np.isnan(subsollon[i]):				# Occasionally we will get a NaN here. This is just for one night -- when Pl was not observed.
    subsollon[i] = 0

  jd[i] = float((words[1])[:-2])                       # Extract JD, and remove comma from string as well

subsollon_rhr = np.round(360 - subsollon,2)     # Need to round this also -- see https://docs.python.org/2/tutorial/floatingpoint.html

##########
# Make a plot to try to match that of Grundy & Fink 1996
##########

rcParams['figure.figsize'] = 18, 10
#if (year_obs == 2014):
#  index_hd = 3		# 0 and 3 are the good HD indices to use
#
#if (year_obs == 2015):
#  index_hd = 1		# Pick a good HD spectrum to use

xlim = (3500, 10000)

if (year_obs == 2014):
  ylim = (0.7, 4)
if (year_obs == 2015):
  ylim = (0.2, 3)

if (year_obs == ''):
  ylim = (0.2, 6)
  
binning = 30
fs      = 24		# Font size

for i in range(len(files_pluto)):
  
#    i = (np.argsort(subsollon_rhr))[j]
    
    is_good = flux_pluto[i] > 0

    plot(wavelength[is_good], smooth_boxcar(flux_pluto[i][is_good] / flux_hd[i][is_good] + i*0.2, binning), color=colors[i])

    plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '' + date_pluto[i] + ', ' + repr(round(subsollon_rhr[i],1)) + '$^\circ$', color=colors[i], fontweight='bold')

plt.title('Pluto/HD, SALT, Compared to Grundy & Fink, binning = ' + repr(binning), fontsize=fs)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Reflectivity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)

file_out = 'spect_ratio_grundyfink_bin' + repr(binning) + '_' + repr(year_obs) + '.png'
plt.savefig(file_out)
print "Wrote: " + file_out
show()



for j in range(len(files_pluto)):
  
    i = (np.argsort(subsollon_rhr))[j]
    print repr(i) + ', subsollon = ' + repr(subsollon_rhr[i]) + ', date_pluto = ' + date_pluto[i]
    
    is_good = flux_pluto[i] > 0

    plot(wavelength[is_good], smooth_boxcar(flux_pluto[i][is_good] / flux_hd[i][is_good] + j*0.2, binning), color=colors[i])

    plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + j*0.2 + 0.75, 
              '' + date_pluto[i] + ', ' + repr(round(subsollon_rhr[i],1)) + '$^\circ$', color=colors[i], fontweight='bold')
show()




##########
# Make a plot of just CH4 band @ 7200
##########

DO_PLOT_7250 = True

if (DO_PLOT_7250):
    xlim = (7000, 7900)
    binning = 30

dy = 0.08

j = 0
lw = 2

#if (year_obs == 2015):
#  ylim = (1.0, 1.8)
#  
#if (year_obs == 2014):
#  ylim = (1.0, 2)

ylim = (1.0, 2)

indices = common_elements(argsort(subsollon), w_is_pluto_good)

for i in indices:
      
    is_good = flux_pluto[i] > 0

    plot(wavelength[is_good], smooth_boxcar(flux_pluto[i][is_good] / flux_hd[i][is_good] + j*dy, binning), 
	  color=colors[mod(j, len(colors))], linewidth=lw)
    plt.text(xlim[1] + 0.02 * (xlim[1]-xlim[0]), 
#   		   flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.05 + 0.75, 
		                                                        j*dy + 1.20, 
		  '  ' + date_pluto[i] + ", " + repr(int(round(subsollon_rhr[i]))) + "$\degree$",
  		  color=colors[j], fontweight='bold', fontsize=18)

    plt.title('Pluto/HD, SALT, CH$_4$ band @ 7250$\AA$, binning = ' + repr(binning) + ', ' + repr(year_obs), fontsize=fs)
    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.ylabel('Reflectivity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)

    j += 1

file_out = 'spect_ratio_7250_bin' + repr(binning) + '_' + repr(year_obs) + '.png'
plt.savefig(file_out)
print "Wrote: " + file_out
show()


##########
# Calculate the sub-solar longitude explicitly. This is as a check. Do it in SPICE.
##########

file_tm  = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"
file_pck = "/Users/throop/gv/dev/kernels/pck00009.tpc"				# Left-handed, so need to do 180- to get RHR
cspice.furnsh(file_tm)
cspice.furnsh(file_pck)

subsollon_spice = np.zeros(len(files_pluto))
name_observer = 'Earth'
name_target   = 'Pluto'
method = 'Near Point'
fixref = 'IAU_PLUTO'
abcorr = 'LT'
r2d    = 360. / 2. / math.pi

for i in range(size(date_pluto)):
    et = cspice.utc2et(date_pluto[i] + ' 20:00:00')				# '18-Aug' means the evening starting 18-Aug
    										# Assume a time of midnight SAST
    spoint = cspice.subsol(method, name_target, et, abcorr, name_observer)
    (radius, lon, lat) = cspice.reclat(spoint)					# Get lon in radians
    if (lon < 0):
      lon += 2. * math.pi
    subsollon_spice[i] = lon*r2d						# Convert to degrees

subsollon_spice_rhr = 360 - subsollon_spice					# This is the final and correct value to plot

##########
# Do some curve fits to the CH4 strength
##########

# We want to fit a linear plus a negative gaussian.

#if (year_obs == 2015):						# 2015: Use full band (** Or for a sanity check, we can use half-band too...
#										# it actually does give similar results)
#alam_0 = 7145      # Wavelength range, in Angstroms
#alam_1 = 7340
#
##if ((year_obs == 2014) or (year_obs == '')):								# 2014: Use abbreviated band (RHS only)
#alam_0_rhs = 7230      # Wavelength range, in Angstroms
#alam_1 = 7340
  
# Set the band wavelength limits. 
# Note that these are slightly different than the limits we use when measuring the area under the curve.
  
alam_gaus     = (7145, 7340)                 # Limits if we are measuring the entire band
alam_gaus_rhs = (7230, 7340)                 # Limits if we are only measuring the narrow ('rhs') side of the CH4 band

m      = 0         # Slope of line
a      = -0.05     # Height of Gaussian
x0     = 7290      # X position of Gaussian
b      = 1.15      # Vertical offset (that is, y-intercept at x=0, *not* at line center)
sigma  = 20        # Width of Gaussian

initial = (m, a, x0, b, sigma)

# First clean up the dataset slightly
# XXX This line is the problem 21-Jul-2015. 'is_good' is different for different spectra. The value here is kind of vestigial, 
# based on what 'i' was used earlier in the code.

rcParams['figure.figsize'] = 10, 5

# Best-fit parameter values
params_fit_gaus_c     = np.zeros((len(files_pluto),3))  # Output parameters for a    constrained gauss + lin fit (3 parameters)
params_fit_gaus_u    = np.zeros((len(files_pluto),5))  # Output parameters for an unconstrained gauss + lin fit (5 parameters)
params_fit_gaus_c_rhs = np.zeros((len(files_pluto),3))  # Output parameters for a    constrained gauss + lin fit (3 parameters). RHS of CH4 band only.
params_fit_gaus_u_rhs = np.zeros((len(files_pluto),5))  # Output parameters for an unconstrained gauss + lin fit (5 parameters). RHS of CH4 band only.

# A best-fit line
line_fit_gaus_c        = np.zeros((len(files_pluto),len(wavelength)))
line_fit_gaus_u        = np.zeros((len(files_pluto),len(wavelength)))
line_fit_gaus_c_rhs    = np.zeros((len(files_pluto),len(wavelength)))
line_fit_gaus_u_rhs    = np.zeros((len(files_pluto),len(wavelength)))

# 
curve_fit_gaus_c        = np.zeros((len(files_pluto),len(wavelength)))
curve_fit_gaus_u        = np.zeros((len(files_pluto),len(wavelength)))
curve_fit_gaus_c_rhs    = np.zeros((len(files_pluto),len(wavelength)))
curve_fit_gaus_u_rhs    = np.zeros((len(files_pluto),len(wavelength)))

# A line that fits the data, used for measuring area under the line.
# This line has different boundaries than the line above (which is used to fit the gaussian)

line_area                   = np.zeros((len(files_pluto), len(wavelength)))
line_area_rhs               = np.zeros((len(files_pluto), len(wavelength)))

# Area under the line
area_line                 = np.zeros(len(files_pluto))
area_line_rhs             = np.zeros(len(files_pluto))

# Perform the fits to the Gaussian + Line. Four fits for each observation.

binning_save = binning
binning = 1    # For the gaussian fits, turn off the binning. It interferes with the NaN's at edges, and also, we want to fit the real data!

for i in range(len(files_pluto)):

    # Full band (wide). These are only going to be valid for the 2015 points -- and not even all of those.
    # Function returns three sets of values... we save this inline. It is perhaps overly compact and very pythonic, I guess.

    (params_fit_gaus_c[i,], success, curve_fit_gaus_c[i,]) = \
        calc_fit_gaus(wavelength, flux_pluto[i] / flux_hd[i], alam_gaus, initial, constrained = True, binning = binning)

    (params_fit_gaus_u[i,], success, curve_fit_gaus_u[i,]) = \
        calc_fit_gaus(wavelength, flux_pluto[i] / flux_hd[i], alam_gaus, initial, constrained = False, binning = binning)
    if (success == False):
        initial = params_fit_gaus_u[i-1,]
        initial = [2.67e-4, -3.24e-2, 7.250e+3, -6.5e-1, 80]
        (params_fit_gaus_u[i,], success, curve_fit_gaus_u[i,]) = \
            calc_fit_gaus(wavelength, flux_pluto[i] / flux_hd[i], alam_gaus, initial, constrained = False, binning = binning)
        
    # Narrow band (RHS only)
    (params_fit_gaus_c_rhs[i,], success, curve_fit_gaus_c_rhs[i,]) = \
        calc_fit_gaus(wavelength, flux_pluto[i] / flux_hd[i], alam_gaus_rhs, initial, constrained = True, binning = binning)

    (params_fit_gaus_u_rhs[i,], success, curve_fit_gaus_u_rhs[i,]) = \
        calc_fit_gaus(wavelength, flux_pluto[i] / flux_hd[i], alam_gaus_rhs, initial, constrained = False, binning = binning)
    
    print 'Finished gaus fit number ' + repr(i) + ', success = ' + repr(success)

    
#    DO_PLOT_INDIVIDUAL_GAUSSIAN = True
#
#    if (DO_PLOT_INDIVIDUAL_GAUSSIAN):
#	plot(data_x, data_y, marker='o', label = 'Data')
#	plot(data_x, out_guess, marker='o', label = 'Initial')
#	plot(data_x, out_line, marker='o', label = 'Initial Linear')
#
#	if (solution):
#	    if (DO_FIT_CONSTRAINED):
#	      out_opt   = gaus_constrained(data_x, popt[0], popt[1], popt[2])		# Gauss( m, a, x0, b sigma): 5 arguments
#	    else:
#   	      out_opt   = gaus(data_x, popt[0], popt[1], popt[2], popt[3], popt[4])	 # 2015: Getting errors here. Why?
#	    out_line2 = popt[0] * data_x + popt[3]
#	    plot(data_x, out_line2, marker='o', label = 'Optimized Linear')
#	    plot(data_x, out_opt, marker='o', label = 'Optimized')
#
#	plt.title("i = " + repr(i) + ", " + date_pluto[i])
#	legend()
#	show()

# Now try a different approach: linear fit and sum under curve
# Left:  7248 (or wherever data start) -- avg over 5 nm starting there
# Left:  7262 (or wherever data start) -- avg over 5 nm starting there
# Right: 7330. This is a pretty clear peak.

#fits_area = np.zeros(len(files_pluto))
#
#if (year_obs == '2015'):				# We have the whole CH4 line to measure -- that's the point
#
#
#if (year_obs == '2014'):				# Only half a spectral line
#  alam_0 = 7262
#  alam_1 = 7330
#

data_x_full = wavelength[is_good]

rcParams['figure.figsize'] = 10, 5

# Now calculate the area under a curve, with endpoints set explicitly
# Turn binning back on for the area curve

binning       = binning_save
alam_area     = (7150, 7330)
alam_area_rhs = (7262, 7330)

for i in range(len(files_pluto)):

    print "Measuring area for curve " + repr(i) + ' with binning ' + repr(binning)    
    (area_line[i], line_area[i,])         = calc_area_curve(wavelength, flux_pluto[i] / flux_hd[i], alam_area,     binning)
    (area_line_rhs[i], line_area_rhs[i,]) = calc_area_curve(wavelength, flux_pluto[i] / flux_hd[i], alam_area_rhs, binning)
    print

##########
# Make a beautiful ganged plot of the individual fit results.
# Sort these in order by longitude, *not* by time.
##########

num_x = 3						# Number of columns to ganged plot
num_y = len(files_pluto)			# Number of rows

#bin_0  = np.min(np.where(data_x_full > alam[0]))
#bin_1  = np.max(np.where(data_x_full < alam[1]))

rcParams['figure.figsize'] = 18,55  # This is the overall size of the whole matrix of plots -- not each individual one
j = 0							# Row number we are currently plotting

fac_expand_range = 1.01

for j in range(len(files_pluto)):

    i = (np.argsort(subsollon_rhr))[j]   # Look up index of the jth element, sorted by longitude
    
    line_y = [1, 1.3] # Y position for vertical lines to draw

# left panel: Area under curves
# Careful here: the the 2014 wide (non-rhs) points are all invalid.
 
    binning = 30
    
    plt.subplot(num_y, num_x, j*num_x+1)		# LHS

    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning))
    plt.ylim((1.05,1.35))

    plot(wavelength, line_area_rhs[i], color='orange', linewidth=1)
    plot(wavelength, line_area[i],     color='orange', linewidth=2) # Wide line for wide range

    plot( [alam_area[0],    alam_area[0]], line_y, color='orange', linewidth=2) # Wide marker for wide interval
    plot( [alam_area[1],    alam_area[1]], line_y, color='orange', linewidth=2)
    plot( [alam_area_rhs[0],alam_area_rhs[0]],  line_y, color='orange', linewidth=1) # Narrow marker for narrow interval
    
    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning), color = 'blue')

    if (j==0):
        plt.title('Area')
        
    plt.xlim((alam_area[0]/fac_expand_range, alam_area[1]*fac_expand_range))

# Middle panel: Gaussian fits, 3-parameter constrained (solid line)

    binning = 1
    
    plt.subplot(num_y, num_x, j*num_x +2)			# Put this plot on RHS, in row j
    
    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning), color = 'blue')
    plt.ylim((1.05,1.35))

    plot(wavelength, curve_fit_gaus_c[i], color = 'lightgreen', linewidth = 2) # Wide line for wide range
    plot(wavelength, curve_fit_gaus_c_rhs[i], color = 'lightgreen', linewidth=1)

    plot( [alam_gaus[0],    alam_gaus[0]], line_y, color='lightgreen', linewidth=2) # Wide marker for wide interval
    plot( [alam_gaus[1],    alam_gaus[1]], line_y, color='lightgreen', linewidth=2)
    plot( [alam_gaus_rhs[0],alam_gaus_rhs[0]],  line_y, color='lightgreen', linewidth=1) # Narrow marker for narrow interval
    
    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning), color = 'blue')

    if (j==0):
        plt.title('Constrained 3-param Gaussian + Linear')
        
    fac_expand_range = 1.01
    plt.xlim((alam_area[0]/fac_expand_range, alam_area[1]*fac_expand_range))

# Right panel: Gaussian Fits, 5-parameter, unconstrained (dashed line)

    binning = 1
    
    plt.subplot(num_y, num_x, j*num_x + 3)			# Put this plot on RHS, in row j
                                                         # plot.subplot(# rows, # columns, Frame # ** starting at 1 **)
    
    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning), color = 'blue')
    plt.ylim((1.05,1.35))
    
    plot(wavelength, curve_fit_gaus_u[i], color = 'red', linewidth = 2) # Wide line for wide range
    plot(wavelength, curve_fit_gaus_u_rhs[i], color = 'red', linewidth = 1)

    plot( [alam_gaus[0],    alam_gaus[0]], line_y, color='red', linewidth=2) # Wide marker for wide interval
    plot( [alam_gaus[1],    alam_gaus[1]], line_y, color='red', linewidth=2)
    plot( [alam_gaus_rhs[0],alam_gaus_rhs[0]],  line_y, color='red', linewidth=1) # Narrow marker for narrow interval
    
    plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i], binning), color = 'blue')

    plt.xlim((alam_area[0]/fac_expand_range, alam_area[1]*fac_expand_range))
    plt.text(7100, 1.3, repr(i) + ', ' + date_pluto[i] + ', ' + repr(round(subsollon_rhr[i],2)) + '${}^\circ$')
    if (j==0):
        plt.title('Unconstrained 5-param Gaussian + Linear')
    
    j += 1
    
plt.show()

# Set appropriate weightings for the three different band depth methods.
# These are entirely arbitrary, just to make the curves fit each other.

weight_fit_area         = 0.2
weight_fit_area_rhs     = 0.95
weight_fit_gaus_u       = 2.5
weight_fit_gaus_c       = 50
weight_fit_gaus_c_rhs   = 40

# Set up flags to indicate which observations are bad, and which are good (e.g., exclude Pluto wide 2014 area measurements, since had no data)
# NB: np.concatenate takes a single argument, which is a list of np arrays. It does not take an arbitrary number of arrays, and cat them.

is_good_area           =  np.concatenate( (np.zeros(13), 1 + np.zeros(7)) ) == 1    # All 2014 bad. All 2015 good. 
                                                                                    # Most 2014 values are NaN anyhow so will not be plotted regardless.  
is_good_area_rhs       = (1 + np.zeros(num_lines)) == 1                      # All of these are good
is_good_gaus_c     = (1 + np.zeros(num_lines)) == 1                      # Green plot, wide band, large dots. All look basically good.
                                                                                    # I probably trust the broad 2014 green more than the rhs 2014 green.
is_good_gaus_c_rhs = (1 + np.zeros(num_lines)) == 1                      # Green plot, narrow band, small dots. All look basically good. 
is_good_gaus_u     = np.concatenate( (np.zeros(18), 1 + np.zeros(2)) ) == 1     # Several just did not fit at all -- not sure why.
is_good_gaus_u_rhs = (1 + np.zeros(num_lines)) == 1                      # All look basically good

is_2014            = np.array(year_pluto) == 2014
is_2015            = np.array(year_pluto) == 2015

# Make a plot of the fits that we have

rcParams['figure.figsize'] = 10,8

w_2014 = np.array(year_pluto) == 2014
w_2015 = np.array(year_pluto) == 2015

plot_elements(subsollon_rhr, area_line_rhs*weight_fit_area_rhs, is_good_area_rhs, linestyle = 'none', marker = 'o', color = 'orange', markersize=5) # narrow line, small symbol
plot_elements(subsollon_rhr, area_line*weight_fit_area,         is_good_area,     linestyle = 'none', marker = 'o', color = 'orange', markersize=10) # wide line, big symbol
plt.ylim((5,0))
plt.title('Area Measurement')
plt.show()

plot(subsollon_rhr, abs(params_fit_gaus_c[:,1]    )*weight_fit_gaus_c,     linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=10)
plot(subsollon_rhr, abs(params_fit_gaus_c_rhs[:,1])*weight_fit_gaus_c_rhs, linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=5) # wide line, big symbol
plt.ylim((0,5))
plt.title('Gaussian fits, 3-parameter constrained')
plt.show()

# Red: 5-param unconstrained fits. There are no wide points here
plot_elements(subsollon_rhr, abs(params_fit_gaus_u[:,1]     * params_fit_gaus_u[:,4])*weight_fit_gaus_u,     is_good_gaus_u    , linestyle = 'none', marker = 'o', color = 'red', markersize=5)
plot_elements(subsollon_rhr, abs(params_fit_gaus_u_rhs[:,1] * params_fit_gaus_u_rhs[:,4])*weight_fit_gaus_u, is_good_gaus_u_rhs, linestyle = 'none', marker = 'o', color = 'red', markersize=10) # wide line, big symbol
plt.ylim((5,0))
plt.title('Gaussian fits, 5-parameter unconstrained')
plt.show()

# Merge all these fits onto one plot -- and plot only the good points

plot_elements(subsollon_rhr, weight_fit_area     * area_line,     is_good_area,     linestyle = 'none', marker = 'o', color = 'orange', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_area_rhs * area_line_rhs, is_good_area_rhs, linestyle = 'none', marker = 'o', color = 'orange', markersize=5)

plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u[:,1]     * params_fit_gaus_u[:,4]    ), is_good_gaus_u,     linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=10)
plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u_rhs[:,1] * params_fit_gaus_u_rhs[:,4]), is_good_gaus_u_rhs, linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=5) # wide line, big symbol

plot_elements(subsollon_rhr, weight_fit_gaus_c     * abs(params_fit_gaus_c[:,1]    ), is_good_gaus_c,     linestyle = 'none', marker = 'o', color = 'red', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_gaus_c_rhs * abs(params_fit_gaus_c_rhs[:,1]), is_good_gaus_c_rhs, linestyle = 'none', marker = 'o', color = 'red', markersize=5) 
plt.ylim((5,0))
plt.xlabel('Sub-Solar Longitude RHR [deg]')
plt.ylabel('CH4 Depth')
plt.title('CH4 7300 Band, all 2014 + 2015 data')

ax = plt.gca()
ax.set_xticks([0, 45, 90, 135, 180, 215, 270, 315, 360])

plt.show()

# Now plot only the 2014 elements

plot_elements(subsollon_rhr, weight_fit_area     * area_line,     np.logical_and(is_good_area, is_2014),     
              linestyle = 'none', marker = 'o', color = 'orange', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_area_rhs * area_line_rhs, np.logical_and(is_good_area_rhs, is_2014), 
              linestyle = 'none', marker = 'o', color = 'orange', markersize=5)

plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u[:,1]     * params_fit_gaus_u[:,4]    ), 
              np.logical_and(is_good_gaus_u, is_2014),    linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=10)
              
plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u_rhs[:,1] * params_fit_gaus_u_rhs[:,4]), 
              np.logical_and(is_good_gaus_u_rhs, is_2014), linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=5) # wide line, big symbol

plot_elements(subsollon_rhr, weight_fit_gaus_c     * abs(params_fit_gaus_c[:,1]    ), 
              np.logical_and(is_good_gaus_c, is_2014),     linestyle = 'none', marker = 'o', color = 'red', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_gaus_c_rhs * abs(params_fit_gaus_c_rhs[:,1]), 
              np.logical_and(is_good_gaus_c_rhs, is_2014), linestyle = 'none', marker = 'o', color = 'red', markersize=5) 
              
plt.ylim((5,0))
plt.xlabel('Sub-Solar Longitude RHR [deg]')
plt.ylabel('CH4 Depth')
plt.title('CH4 7300 Band, 2014 only')
plt.show()

# Now plot only the 2015 elements

plot_elements(subsollon_rhr, weight_fit_area     * area_line,     np.logical_and(is_good_area, is_2015),     
              linestyle = 'none', marker = 'o', color = 'orange', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_area_rhs * area_line_rhs, np.logical_and(is_good_area_rhs, is_2015), 
              linestyle = 'none', marker = 'o', color = 'orange', markersize=5)

plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u[:,1]     * params_fit_gaus_u[:,4]    ), 
              np.logical_and(is_good_gaus_u, is_2015),    linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=10)
              
plot_elements(subsollon_rhr, weight_fit_gaus_u * abs(params_fit_gaus_u_rhs[:,1] * params_fit_gaus_u_rhs[:,4]), 
              np.logical_and(is_good_gaus_u_rhs, is_2015), linestyle = 'none', marker = 'o', color = 'lightgreen', markersize=5) # wide line, big symbol

plot_elements(subsollon_rhr, weight_fit_gaus_c     * abs(params_fit_gaus_c[:,1]    ), 
              np.logical_and(is_good_gaus_c, is_2015),     linestyle = 'none', marker = 'o', color = 'red', markersize=10) # wide line, big symbol
plot_elements(subsollon_rhr, weight_fit_gaus_c_rhs * abs(params_fit_gaus_c_rhs[:,1]), 
              np.logical_and(is_good_gaus_c_rhs, is_2015), linestyle = 'none', marker = 'o', color = 'red', markersize=5) 
              
plt.ylim((5,0))
plt.xlabel('Sub-Solar Longitude RHR [deg]')
plt.ylabel('CH4 Depth')
plt.title('CH4 7300 Band, 2015 only')
plt.show()

# Now average all the elements, night-by-night

strength = []
stdev    = []

for i in range(num_lines):
  strength_i = np.array([\
  np.nan_to_num(weight_fit_area * area_line * is_good_area)[i], \
  np.nan_to_num(weight_fit_area_rhs * area_line_rhs * is_good_area_rhs)[i], \
  np.nan_to_num(weight_fit_gaus_u * abs(params_fit_gaus_u[:,1]     * params_fit_gaus_u[:,4]    ) * is_good_gaus_u)[i], \
  np.nan_to_num(weight_fit_gaus_u * abs(params_fit_gaus_u_rhs[:,1]     * params_fit_gaus_u_rhs[:,4]    ) * is_good_gaus_u_rhs)[i], \
  np.nan_to_num(weight_fit_gaus_c * abs(params_fit_gaus_c[:,1]         ) * is_good_gaus_c)[i], \
  np.nan_to_num(weight_fit_gaus_c * abs(params_fit_gaus_c_rhs[:,1]     ) * is_good_gaus_c_rhs)[i]])
  
  elements_valid_i = (strength_i != 0)
  stdev_i = np.std(strength_i[elements_valid_i])
  strength.append(np.sum(strength_i) / np.sum((strength_i !=0)*1))
  stdev.append(stdev_i)

# Calculate the final band strength, and the error bar (which is from the variation in our measurements of each point)

strength = np.array(strength)
stdev = np.array(stdev)

##########
# Make a 'final, everything' plot of band strength vs. longitude
##########

plot(subsollon_rhr, strength, marker = 'o', linestyle='none', markersize=6)
plt.errorbar(subsollon_rhr, strength, yerr=stdev, linestyle='none', linewidth=1.5)

plt.ylim((5,0))
ax = plt.gca()
ax.set_xticks([0, 45, 90, 135, 180, 215, 270, 315, 360])

plt.xlabel('Subsolar Longitude [RHR]', fontsize=15)
plt.ylabel('Strength', fontsize=15)
plt.xlim((0,360))
plt.title('CH$_4$ band strength, 730 nm', fontsize=18)

plt.show()

##########
# Plot center position of CH4 band vs. longitude
# For this use only the unconstrained Gaussian fits
##########

plot(subsollon, params_fit_gaus_u[:,2], linestyle='none', marker = 'o', color = 'lightgreen', markersize=10)
plot(subsollon, params_fit_gaus_u_rhs[:,2], linestyle='none', marker = 'o', color = 'lightgreen', markersize=5)

plt.xlabel('Subsolar longitude [RHR]', fontsize=20)
plt.ylabel('Band position [$\AA$]', fontsize=20)
plt.ylim((7240, 7310))
ax = plt.gca()
ax.set_xticks([0, 45, 90, 135, 180, 215, 270, 315, 360])
plt.xlim((0,360))
plt.show()

##########
# Plot center position of CH4 band vs. longitude
# For this use only the unconstrained Gaussian fits
##########

#plot(subsollon, params_fit_gaus_u[:,4], linestyle='none', marker = 'o', color = 'lightgreen', markersize=10)

plot(subsollon, params_fit_gaus_u_rhs[:,4], linestyle='none', marker = 'o', color = 'lightgreen', markersize=5)
plt.xlabel('Subsolar longitude [RHR]', fontsize=20)
plt.ylabel('Band Width [$\AA$]', fontsize=20)
plt.ylim((5,28))
ax = plt.gca()
ax.set_xticks([0, 45, 90, 135, 180, 215, 270, 315, 360])
plt.xlim((0,360))
plt.show()

##########
# Plot an overall spectrum for 2014 and for 2015
# We have already excluded bad spectra -- all of these should be good.
##########

rcParams['figure.figsize'] = 16,6

binning = 2

flux_pluto_arr = np.array(flux_pluto)
flux_hd_arr    = np.array(flux_hd)

# Convert all zero to nan. All invalid numbers are supposed to be nan already, but some 0's crept in.
# We want nan since np.nanmean() will skip those elements, and mean over the rest.

flux_pluto_arr[where(flux_pluto_arr == 0)] = nan
flux_hd_arr[where(flux_hd_arr == 0)] = nan

# There are a lot of different ways to combine using sum, mean, nanmean, etc.
# I think this is the best one
# Ugh... all good, except it doens't really work. Ch4 band is messed up.

iof_2015 =  np.nanmean( flux_pluto_arr[is_2015,:] / flux_hd_arr[is_2015,:],0)
iof_2014 =  np.nanmean( flux_pluto_arr[is_2014,:] / flux_hd_arr[is_2014,:],0)

flux_pluto_all = np.nanmean

# Try two different ways of combining the data. They have the same shape and noise, but are more different than I would think.
# The only times these should differ are for regions that have NaN in the HD and data in the Pluto, or v/v. And this is not very much!

iof_all_1  =  np.nanmean( flux_pluto_arr / flux_hd_arr, 0) # This is probably technically more accurate -- we average the I/F's.
iof_all_2  =  np.nanmean( flux_pluto_arr,0) / np.nanmean(flux_hd_arr, 0) # This will improve SNR, 

iof_all = iof_all_1 # Pick with one and go with it.

# Do one minor hand-edit to the spectrum

iof_all[7260:7269] = nan # Get rid of one spike.

rcParams['figure.figsize'] = 16,6
plot(wavelength, smooth_boxcar(iof_all,30))
plt.xlim((3400,9400))
plt.ylim((0.35, 1.45))
plt.xlabel('Wavelength [$\AA$]', fontsize=20)
plt.title('I/F SALT, All Spectra 2014-2015', fontsize=20)
plt.ylabel('I/F [arbitrary]', fontsize=20)

#plot(wavelength, smooth_boxcar(spect_2015*1.4, binning), color = 'red')
plot(wavelength, smooth_boxcar(flux_2014,     binning), color = 'orange')
plot(wavelength, smooth_boxcar(flux_2015_almost*1.75,     binning), color = 'lightgreen')

# For debugging: quickly plot all 2015 spectra

quit

for i in nditer(where(is_2015)): # nditer is a 'generating function' that loops over every value of where().
    plot(wavelength, smooth_boxcar(flux_pluto_arr[i] / flux_hd_arr[i] + i*0.2,30))
plot(wavelength, iof_all_1+2, color = 'black')    
plt.title('Pluto 2015')    
plt.show()

for i in nditer(where(is_2015)): # nditer is a 'generating function' that loops over every value of where().
    plot(wavelength, flux_hd[i] + i*0.2)
plt.title('HD 2015')    
plt.show()

rcParams['figure.figsize'] = 16,12

for i in range(num_lines):
#    plot(wavelength, flux_pluto_arr[i,:] / flux_hd_arr[i,:] + i * 0.2)
   plot(iof_arr[i,:] + i * 0.2)

plot(iof_all - 0.2)
plt.xlim((8400,8500))
    
# Get a total I/F, merging all spectra and dividing by all HD

flux_pluto_denan = np.copy(flux_pluto)
flux_hd_copy = np.copy(flux_hd)
flux_pluto_denan = np.nan_to_num(f)

val = nan
iof_all[3145] = val
iof_all[3287:3306] = val
iof_all[3424:3425+1] = val
iof_all[5405:5410+1] = val
iof_all[5518:5522+1] = val
iof_all[5676:5678]   = val
iof_all[7260:7268]   = val
iof_all[8339:8344+1] = val
iof_all[10364:10373] = val
iof_all[12288:12299] = val

#plot(wavelength, iof_all, markersize = 1, linestyle = 'none', marker='.')n
plot(iof_all, markersize = 1, linestyle = 'none', marker='.')
plot(iof_all, markersize = 1, linestyle = 'none', marker='.')
#plt.xlim((3000,6000))



#plt.xlim((3000,4000))


flux_
# Remove a few spikes

# Add the two spectra for a final spectrum

plot(wavelength, smooth_boxcar(flux_2014 + flux_2015, 1))

# Plot the 620 nm CH4 feature

rcParams['figure.figsize'] = 6,6
plot(wavelength, smooth_boxcar(flux_2014 + flux_2015_almost, 1))
plt.xlim((6000,6500))
plt.ylim((15,25))

p0 =  (m, a, x0, b, sigma)
m = 1/500.   # Rough slope
b = 1000    # Rough intercept
x0 = 6200
sigma = 8
a = -0.1      # Band depth. I assume sign is right?
initial = (m, a, x0, b, sigma)
alam_6200 = ((6000,6500))

(fits,success,data_fit_out) = calc_fit_gaus(wavelength, flux_2014_2015, alam_6200, initial, constrained=False, binning=1)

plot(wavelength, flux_2014_2015, color='lightgreen')
plot(wavelength, data_fit_out, color='red')
plt.xlim((6000,6500))
plt.ylim((20,24))
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('I/F')
plt.title('Pluto, 620 nm band')



plt.xlim((3400,9500))
plt.ylim((4, 19))
plt.xlabel('Wavelength [$\AA$]', fontsize=20)
plt.ylabel('Intensity', fontsize=20)
plt.title('Pluto, 2014 and 2015')
plt.show()

# Plot 
quit

# recompute the ax.dataLim
#    ax.relim()
## update ax.viewLim using the new dataLim
#    ax.autoscale_view(True,True,True)
#    plt.draw()
#    plt.show()


#    plt.title('i = ' + repr(i) + ', subsollon = ' + repr(subsollon[i]) + ', date = ' + date_pluto[i])
#    plt.xlabel('Wavelength')
#    plt.ylabel('Intensity')
#    show()
    


# Calculate a slope of the continuum in the red, and in the blue.
# I have just chosen a few wavelength ranges, essentially at random, because they looked pretty clean
#
# REMOVING ALL SPECTRAL SLOPE THINGS SINCE I THINK THEY'RE UNRELIABLE
#
#   4000-4300
#   5800-6300
#
#alam_slope1 = [4000, 4300]
#alam_slope2 = [5950, 6800]
#alam_slope3 = [4700, 5300]
#
#bin_slope1 = [np.min(np.where(data_x_full > alam_slope1[0])), 
#              np.min(np.where(data_x_full > alam_slope1[1]))] 
#
#bin_slope2 = [np.min(np.where(data_x_full > alam_slope2[0])), 
#              np.min(np.where(data_x_full > alam_slope2[1]))] 
#
#bin_slope3 = [np.min(np.where(data_x_full > alam_slope3[0])), 
#              np.min(np.where(data_x_full > alam_slope3[1]))] 
#
#bin_1  = np.max(np.where(data_x_full < alam_1))
#
#slope1     = np.zeros(len(files_pluto))
#slope2     = np.zeros(len(files_pluto))
#slope3     = np.zeros(len(files_pluto))
#intercept1 = np.zeros(len(files_pluto))
#intercept2 = np.zeros(len(files_pluto))
#intercept3 = np.zeros(len(files_pluto))
#
#for i in range(len(files_pluto)):
#
#    data_y_full = (flux_pluto[i] / flux_hd[3])[is_good]
#    data_x_full = wavelength[is_good]
#
#    data_x = data_x_full[bin_slope1[0]:bin_slope1[1]]
#    data_y = data_y_full[bin_slope1[0]:bin_slope1[1]]
#
#    mdx = np.mean(data_x)			# Mean of all x values
#    mdy = np.mean(data_y)			# Mean of all y values
#
## Use the formula for linear least-squares fit to get the slope
#
#    slope1[i] = np.sum( (data_x - mdx) * (data_y - mdy) ) / np.sum( (data_x - mdx)**2)
#    intercept1[i] = mdy - slope1[i] * mdx
#
## Same for slope2
#
#    data_x = data_x_full[bin_slope2[0]:bin_slope2[1]]
#    data_y = data_y_full[bin_slope2[0]:bin_slope2[1]]
#
#    mdx = np.mean(data_x)
#    mdy = np.mean(data_y)
#
#    slope2[i] = np.sum( (data_x - mdx) * (data_y - mdy) ) / np.sum( (data_x - mdx)**2)
#    intercept2[i] = mdy - slope2[i] * mdx
#
## Same for slope3
#
#    data_x = data_x_full[bin_slope3[0]:bin_slope3[1]]
#    data_y = data_y_full[bin_slope3[0]:bin_slope3[1]]
#
#    mdx = np.mean(data_x)
#    mdy = np.mean(data_y)
#
#    slope3[i] = np.sum( (data_x - mdx) * (data_y - mdy) ) / np.sum( (data_x - mdx)**2)
#    intercept3[i] = mdy - slope3[i] * mdx


##########
# Make Plot of sub-solar longitude vs. slopes
##########

rcParams['figure.figsize'] = 12,7 

w = is_pluto_good		# Which element to plot.

if (year_obs == 2014):
  ylim = (1,3)
if (year_obs == 2015):
  ylim = (0.5,5.5)

fac = 10000			# Multiply the data by this scaling factor

subsollon_0 = np.copy(subsollon_rhr)	# This particular plot shows a trend centered on longitude 0.  Define a new longitude system like that
subsollon_0 -= 360 * (subsollon_0 > 180)

plot(subsollon_0[w], slope1[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope1) + "$\AA$", color = 'blue')
plot(subsollon_0[w], slope3[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope3) + "$\AA$", color = 'green')
plot(subsollon_0[w], slope2[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope2) + "$\AA$", color = 'red')

#   for i in range(len(subsollon)):
#     if is_pluto_good[i]:
#       plt.text(subsollon_0[i], 5, '  ' + date_pluto[i][5:], fontsize = fs)
  
plt.ylim(ylim)
plt.title('Pluto Spectral Slope vs. Longitude, SALT ' + repr(year_obs), fontsize = fs * 2)
plt.xlabel('Pluto Sub-Solar Longitude', fontsize = fs * 2)
plt.ylabel('Slope $(10^4 \AA)^{-1}$', fontsize = fs*2)
plt.xticks([-180, -90, 0, 90, 180])
legend()
show()

##########
# Make Plot of sub-solar longitude vs. slopes
##########

rcParams['figure.figsize'] = 12,7 

w = is_pluto_good		# Which element to plot.

fac = 10000			# Multiply the data by this scaling factor

subsollon_0 = np.copy(subsollon_rhr)	# This particular plot shows a trend centered on longitude 0.  Define a new longitude system like that
subsollon_0 -= 360 * (subsollon_0 > 180)

#   plot(subsollon_0[w], slope1[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope1) + "$\AA$", color = 'blue')
#   plot(subsollon_0[w], slope3[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope3) + "$\AA$", color = 'green')
plot(subsollon_0[w], slope2[w]*fac, marker = 'o', markersize=ms, linestyle='none', label = range2str(alam_slope2) + "$\AA$", color = 'red')

#   for i in range(len(subsollon)):
#     if is_pluto_good[i]:
#       plt.text(subsollon_0[i], 5, '  ' + date_pluto[i][5:], fontsize = fs)
  
plt.ylim((1,2))
plt.title('Pluto Spectral Slope vs. Longitude, SALT ' + repr(year_obs), fontsize = fs * 2)
plt.xlabel('Pluto Sub-Solar Longitude', fontsize = fs * 2)
plt.ylabel('Slope $(10^4 \AA)^{-1}$', fontsize = fs*2)
plt.xticks([-180, -90, 0, 90, 180])
legend()
show()

quit

##########
# Make a plot to try to match that of Grundy & Fink 1996 -- but over the whole wavelength range
##########
              
# Align on some feature near the middle of wavelength range. Looks like a good one at ~ 6870

align_wavelengths(wavelength, flux_pluto, (6850, 6890), 6865, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (6850, 6890), 6865, smoothing=10)
  
xlim = (3200, 10000)

ylim = (0.3, 4)
binning = 30

for i in range(len(files_pluto)):

      is_good = flux_pluto[i] > 0

      plot(wavelength[is_good], smooth_boxcar(flux_pluto[i][is_good] / flux_hd[i][is_good] + i*0.2, binning), color=colors[i])

      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])
              
plt.title('Pluto/HD, SALT, binning = ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)

file_out = 'spect_ratio_full_{0}-{1}_bin{2}.png'.format(xlim[0], xlim[1], binning)
plt.savefig(file_out)
print "Wrote: " + file_out

show()

##########
# Make a plot of 5500 - 6500 to try to see O2 bands of Calvin & Spencer
# Should see bands at 5600-5700 and 6200 - 6300.
# SHould see a bit peak up at 5750. I don't see it...
##########
              
# Align on some feature near the middle of wavelength range. Looks like a good one at ~ 6870

align_wavelengths(wavelength, flux_pluto, (6850, 6890), 6865, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (6850, 6890), 6865, smoothing=10)
  
xlim = (5500, 6500)

ylim = (0.3, 4)
binning = 1

for i in range(len(files_pluto)):
      plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i] + i*0.2, binning), color=colors[i])
      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])
              
plt.title('Pluto/HD, SALT, binning = ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)

file_out = 'spect_ratio_full_{0}-{1}_bin{2}.png'.format(xlim[0], xlim[1], binning)
plt.savefig(file_out)
print "Wrote: " + file_out

show()

# Make plots of region around 0.89 um band
# Loop over this, trying each HD observation, to find the best one for each Pluto spectrum
# Make 12 plots, each with 5 lines

#rcParams['figure.figsize'] = 20, 6
#
#xlim = (8300, 9300)
#index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy. Does it really go to 14?
#for i in range(len(files_pluto)):
#  if i in index_good_spectra_pluto:
#      for j in range(len(files_hd)):
#        plot(wavelength, flux_pluto[i] / flux_hd[j] + j*0.3, color=colors[j], 
#             label = 'Pluto #' + repr(i) + ' / HD #' + repr(j))
#
#      plt.title('Pluto #' + repr(i) + ' / HD #n', fontsize=24)
#      plt.xlim(xlim)
#      plt.ylabel('Intensity [arbitrary]', fontsize=24)
#      plt.xlabel('Wavelength [$\AA$]', fontsize=24)
#      plt.ylim((1.0,3))
#      legend()
#      show()

# Make plot with all Pluto / HD, scaled to match that in Grundy paper

xlim = (5000, 9300)
binning = 5
# rcParams['figure.figsize'] = 12, 6

index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, smooth_boxcar(flux_pluto[i] / flux_hd[i] + i*0.2, binning), color=colors[i])
      plt.text(amax(wavelength), 
               np.median(flux_pluto[i][-200:-50]) + i*0.2, 
               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT ' + repr(year_obs), fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim((0.7, 4))

# Now make plot of the ratios of the HD stars night-to-night
# This shows us the atmospheric variability.
# Concl: it looks like a few percent

##########
# HD stars, mutual ratios, narrow wavelength range
##########

# xlim = (8300, 9300)
# for i in range(len(files_hd)):
#       for j in range(len(files_hd)):
#         plot(wavelength, flux_hd[i] / flux_hd[j] + j*1, color=colors[j], 
#              label = 'HD #' + repr(i) + ' / HD #' + repr(j))
# 
#       plt.title('HD #' + repr(i) + ' / HD #n', fontsize=24)
#       plt.xlim(xlim)
#       plt.ylabel('Intensity [arbitrary]', fontsize=24)
#       plt.xlabel('Wavelength [$\AA$]', fontsize=24)
#       plt.ylim((0.5,5.5))
#       legend()
#       show()

# Now make plot of the ratios of the HD stars night-to-night
# This shows us the atmospheric variability.
# ** Same as above, but different wavelength range

##########
# HD stars, mutual ratios, full wavelength range
##########

# xlim = (3300, 9300)
# for i in range(len(files_hd)):
#       for j in range(len(files_hd)):
#         plot(wavelength, flux_hd[i] / flux_hd[j] + j*1, color=colors[j], 
#              label = 'HD #' + repr(i) + ' / HD #' + repr(j))
# 
#       plt.title('HD #' + repr(i) + ' / HD #n', fontsize=24)
#       plt.xlim(xlim)
#       plt.ylabel('Intensity [arbitrary]', fontsize=24)
#       plt.xlabel('Wavelength [$\AA$]', fontsize=24)
#       plt.ylim((0.5,5.5))
#       legend()
#       show()

# Make a plot of airmasses vs. date, for both HD and Pluto

plot(jd, subsollon, color= 'blue', label = 'Pluto', linestyle = 'none', marker = 'o')
#plot(jd[is_hd_good],    airmass[is_hd_dood],    color='red',   label = 'HD',    linestyle = 'none', marker = 'o')

plt.ylabel('Subsollon')
plt.xlabel('JD')
plt.legend()
show()

# Measure the slope of the blue spectrum
# Do this by fitting a line to the 

# Now experiment with binning

plot(wavelength, smooth_boxcar(flux_pluto[2] / flux_hd[4], 2))
plt.ylim((0.4, 1.4))
plt.xlim((8300,9300))
show()

rcParams['figure.figsize'] = 10, 28

plt.figure(1)
plt.subplot(12,2,1)		# Can abbreviate (2,2,3) as (223). nrows, ncolums, plot_#. Latter is counted from UL, to LR
plt.plot([4,5])

plt.subplot(12,2,2)
plt.plot([4,7])

plt.subplot(12,2,10)
plt.plot([4,7])

plt.subplot(12,2,20)
plt.plot([4,7])


show()

#########
# Function for area_line()
#########


def area_line ( data_x_full, data_y_full, x_0, x_1, x_0_plot, x_1_plot, title='', do_plot=False ):
    "Determines the area between a curve (defined by data_) and a line (defined by xy)"

    rcParams['figure.figsize'] = 4, 02

    bin_0  = np.min(np.where(data_x_full > x_0))
    bin_1  = np.max(np.where(data_x_full < x_1))

    bin_0_plot  = np.min(np.where(data_x_full > x_0_plot))
    bin_1_plot  = np.max(np.where(data_x_full < x_1_plot))

#       plt.subplot(num_y, num_x, j*num_x)			# Put this plot on RHS, in row j

    data_y = data_y_full[bin_0:bin_1]
    data_x = data_x_full[bin_0:bin_1]

    data_y_plot = data_y_full[bin_0_plot:bin_1_plot]
    data_x_plot = data_x_full[bin_0_plot:bin_1_plot]

# Extract just the relevant portion

    dy = data_y[-1] - data_y[0]
    dx = data_x[-1] - data_x[0]

    m = dy / dx		      # Slope
    b = data_y[0] - m * data_x[0] # Y-intercept

    line = m * data_x + b				# Line

    area = sum(line - data_y) 		# Sum the total area. Line - data because data is below. Don't use ABS().

    if (do_plot):
	plot(data_x, line)					# Plot the line 
	plot(data_x, data_y)				# Plot the data
	plot(data_x_plot, data_y_plot)			# Plot the data
#   	plt.title(title + ', area = ' + repr(area))
#   	plt.title(title)
#           plt.txt(7050, (np.max(data_y) - np.min(data_y))/2, "AA")
	show()
    return area

##########
# Some testing for plot_area
##########

##########
# First, measure just half the band (that is, the RHS of it)
##########

rcParams['figure.figsize'] = 12,7

alam_0      = 7250      # Wavelength range, in Angstroms
alam_1      = 7330
alam_0_plot = 7100
alam_1_plot = 7500
index_hd    = 0		# 1/2/3 are consistent. 4 is way off. 0 is consistent but requires +5 offset.
fit_area    = []
fs 	    = 24

for i in range(size(subsollon)):		# 'Returns indices that would sort an array'

    data_x_full = wavelength[is_good]		# data_{xy}_good are where there is non-zero flux

    data_y_full = smooth_boxcar(flux_pluto[i][is_good] / flux_hd[index_hd][is_good] + i*0.05, binning)

    title = 'i = ' + repr(i) + ', subsollon_rhr = ' + repr(subsollon_rhr[i])

    fit_area.append( area_line ( data_x_full, data_y_full, alam_0, alam_1, alam_0_plot, alam_1_plot, title=title, do_plot=False ) )

#   plot(subsollon_rhr, fit_area, linestyle='none', marker = 'o')
#   plt.ylim(5, -1)

fit_area_narrow = np.array(fit_area)

if (year_obs == 2014):				# We could make this plot for either 2014 or 2015 (that is, full CH4 or narrow CH4)

# Delete three data points. I think these are from bad spectra. I deleted them from my NH SALT talk, so I'll delete them here too.

    fit_area[12] = -1
    fit_area[6]  = -1
    fit_area[5]  = -1

    plot(subsollon_rhr, fit_area, linestyle='none', marker = 'o', label = 'Narrow', ms = ms, color = 'blue')

    plt.title('SALT Pluto 7300$\AA$ CH$_4$ Depth vs Longitude, HD[' + repr(index_hd) + '], ' + repr(year_obs), fontsize=fs)
    plt.xlabel('Subsolar Longitude [deg, RHR]', fontsize=fs)
    plt.ylabel('Band Depth', fontsize=fs)
    #   plt.title('Pluto 7300$\AA$ CH$_4$ Depth vs. Longitude, SALT ' + repr(year_obs), fontsize = fs*2)
    plt.xticks([0, 90, 180, 270, 360])
    plt.ylim(6, 0)
    plt.xlim(-10, 370)

    file_out = 'depth_7300_v_longitude_area_' + repr(year_obs) + '.png'
    plt.savefig(file_out)
    print "Wrote: " + file_out
    show()
  
##########
# Calculate and make a plot of the CH4 strength vs. longitude, using area method, for wide CH4 band
# To run this, run the narrow method first (above), and then run this -- otherwise will get broadcasting error.
##########

if (year_obs == 2015):		# This only makes sense to do for 2015, where we can use the full band

    index_hd    = 2		# 1/2/3 are consistent. 4 is way off. 0 is consistent but requires +5 offset.
    alam_0      = 7100      # Wavelength range, in Angstroms
    #   alam_1 = 7330	# Nominal RHS of band
    alam_1      = 7375 		# Very wide -- right-most possible edge, a bit beyond the edge
    alam_0_plot = 7000
    alam_1_plot = 7600

    fit_area = []

    binning = 30 

    for i in range(size(subsollon)):		# 'Returns indices that would sort an array'

	data_x_full = wavelength[is_good]		# data_{xy}_good are where there is non-zero flux

	data_y_full = smooth_boxcar(flux_pluto[i][is_good] / flux_hd[index_hd][is_good] + i*0.05, binning)

	title = 'i = ' + repr(i) + ', subsollon_rhr = ' + repr(subsollon_rhr[i])
	title = repr(int(subsollon_rhr[i])) + "$\degrees$"

	fit_area.append( area_line ( data_x_full, data_y_full, alam_0, alam_1, alam_0_plot, alam_1_plot, title=title, do_plot=True ) )

    #   plot(subsollon_rhr, fit_area, linestyle='none', marker = 'o')
    #   plt.ylim(10, -1)
    show()

    fit_area_wide = np.array(fit_area)

    ms = 10 # marker size
    fs = 24 # font size

## CHANGE ONE DATA POINT
    fit_area_wide[0] = fit_area_narrow[0]*3 			# This is to correct one data point with a missing error bar
    fit_area_wide[0] = 3 			# This is to correct one data point with a missing error bar

    fit_area_mean = ( fit_area_wide+5 + fit_area_narrow*3 ) / 2

    #   plot(subsollon_rhr, fit_area_wide+5, linestyle='none', marker = 'o', label = 'Wide', ms = ms, color = 'red')
    #   plot(subsollon_rhr, fit_area_narrow*3, linestyle='none', marker = 'o', label = 'Narrow', ms = ms, color = 'blue')
    plot(subsollon_rhr, fit_area_mean, linestyle='none', marker = 'o', label = 'Narrow', ms = ms, color = 'red')


    yerr = np.abs(fit_area_wide - fit_area_narrow*3) / 2


    plt.title('SALT Pluto 7300$\AA$ CH$_4$ Depth vs Longitude, HD[' + repr(index_hd) + '], ' + repr(year_obs), fontsize=fs)
    plt.xlabel('Subsolar Longitude [deg, RHR]', fontsize=fs)
    plt.ylabel('Band Depth', fontsize=fs)
    #   plt.title('Pluto 7300$\AA$ CH$_4$ Depth vs. Longitude, SALT ' + repr(year_obs), fontsize = fs*2)
    plt.xticks([0, 90, 180, 270, 360])
    plt.ylim(13, -1)
    plt.xlim(-10, 370)
    plt.errorbar(subsollon_rhr, fit_area_mean, xerr=0, yerr=yerr, ls = 'none', color = 'red', linewidth=3)


    file_out = 'depth_7300_v_longitude_area_' + repr(year_obs) + '.png'
    plt.savefig(file_out)
    print "Wrote: " + file_out
    show()


#   fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)
#   ax = axs[0,0]
#   ax.errorbar(x, y, yerr=yerr, fmt='o')
#   ax.set_title('Vert. symmetric')

