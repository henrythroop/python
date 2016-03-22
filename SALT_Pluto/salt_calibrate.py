# These are batch commands to calibrate SALT data using PyRAF / PySALT

# Program takes the raw FITS files, and applies the SALT pipeline to them.
# Output is .txt file for each spectrum, one per file.
# It does not combine these or make (non-trivial) plots.
#
# HBT 15-Aug-2014

# Set the prompt color. Can't figure out how to put magic command into a config file in .ipython...

# %config PromptManager.in_template = r'{color.LightGreen}In [{color.LightGreen}{count}{color.LightGreen}]: '

from pyraf import iraf
from iraf import pysalt
from iraf import saltspec
from iraf import specextract
from iraf import specprepare
from iraf import specidentify
from iraf import specrectify
from iraf import specsky
from iraf import saltred

import sys # for sys.sdout.write()

import pdb
import subprocess

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
import astropy
from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 
from astropy.convolution import convolve, Box1DKernel

import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *

from matplotlib.path import Path
import matplotlib.patches as patches

##########
def shorten_filename(file):
##########
    out = file.split('/')[-1]  # Get just the filename itself
    out = out.split('_')[0]  # Get just the filename itself
    return out

##########
def get_section_spect(file, dy_spect):		# Return a string [100:120] stating which region of the image to extract
##########
      row_max = get_row_max(file)

      section = '[' + repr(row_max + -dy_spect/2) + ':' + repr(row_max + dy_spect/2) + ']'

      return section

##########
def get_section_sky(file, y0_sky, dy_sky):		# returns a string, like "[1:100]"
##########

    hdulist = fits.open(file)
    image = hdulist['SCI'].data

    row_max = get_row_max(file)
    section_sky = "[{0:d}:{1:d}]".format(int(row_max + y0_sky), int(row_max + y0_sky + dy_sky))
    hdulist.close()

    return section_sky

##########
def get_year_obs(file):
##########
  year_obs = 0
  if (file.find('P2015') > 1): year_obs = 2015
  if (file.find('P2014') > 1): year_obs = 2014

  return year_obs

##########
def get_row_max(file):		# Look at a spectrum and return what row the brightest signal is at
##########
  
       if (os.path.isfile(file) == 0):		# Return 0 if file doesn't exist
         return 0

       hdulist = fits.open(file)
       image = hdulist['SCI'].data

       year_obs = get_year_obs(file)

#   y0_sky[ np.char.find(files, '06280032') > 0 ] = 40
   
   # Look up the brightest row. Or something close to that, using a more intelligent algorithm.
   # ** For observations on 20140816, there are two objects in the field. For these, we must limit the search
   #    width to +- 90 pixels; otherwise, the star will be chosen, not Pluto.

# 2015 data is 684 rows=spatial x 3171 columns = wavelength

# Handle Year 2014

       if (year_obs == 2014):
         row_center = 1022               # Row that the center of the spectrum is often on for RSS
         row_min    = row_center - 200   # Search only near the center.
         row_max    = row_center + 200
   
         if (file.find('20140816') > 1): # Special case here since bright star in nearby field
             row_min = row_center - 100
             row_max = row_center + 100

# Handle Year 2015 (fewer rows)

       if (year_obs == 2015):
         row_center = 340
         row_min    = row_center - 50   # Search only near the center.
         row_max    = row_center + 50
           
       rowsum = image.sum(1)
       
   # Set it to zero outside the center, so to limit the search to bright lines near the center
   
       rowsum[0:row_min-1] = 0
       rowsum[row_max:] = 0
       
   # Smooth it
   
       height_smooth = 11 # Must be odd
       rowsum = convolve(rowsum, Box1DKernel(height_smooth)) # Convolving it shifts it because output array is longer than input
       rowsum = roll(rowsum, -(height_smooth-1)/2)           # Shift it back...
   #    rowsum = rowsum - np.median(rowsum) # subtract off a background
     
       row_max = where(rowsum == amax(rowsum))[0][0]
   
       hdulist.close()

       return row_max

##########
# Print status line for this file
##########

def print_status_line_file(j, file, object, exptime, specid, cr, rect, sky, ext):

    print "{0:3d}. {1} {2:10} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d} {8:5d} {9:5d}.".format(
	j, file, object, exptime, specid,
	cr, rect, sky, ext, j)

##########
# Make a plot of a spectrum
##########

def plot_spectrum(file, object, date):
    print "loading spectrum from file " + file
    d = loadtxt(file)
    wavelength = d[:,0]
    flux = d[:,1]
    dflux = d[:,2]
    yrange = ((np.mean(flux) - 2*np.std(flux), np.mean(flux) + 2*np.std(flux)))

    rcParams['figure.figsize'] = 10, 6 # Make plot normal size
#       colors = np.array(['grey', 'brown', 'orange', 'yellow', 'blue', 'red', 'pink', 'green'])
    plot(wavelength,flux, color='red', ls='-')
    plt.title(shorten_filename(file) + ' ' + object + ' ' + date)
    plt.ylim(yrange)
    plt.show()

#       for i in range(len(fluxes)):
#   	if (object == 'Pluto'):
#   	    f = np.array(flux).astype(float)
#   	    f = f / np.amax(f)  # normalize it

#       plt.axis([3000,9000,0,100])
#       plt.xlim((3000,9000))
        
##########
# Make a plot of the image rectified
##########

def plot_rectify(file, object, date):

    hdulist = fits.open(file)
    image = hdulist['SCI'].data
    dy_spect_plot = 400		# How many rows vertically is the spectrum? ** Changes 2014 vs. 2015
    ny = image.shape[0]	# Look up the central row of the spectrum, and take +- dy/2 rows from that
    plt.imshow(log(image))
    plt.ylim((ny/2 - dy_spect_plot/2, ny/2 + dy_spect_plot/2))

    plt.title(file + ", " + object + ", " + ' RECT')
    plt.show()
    hdulist.close()

##########
# Make a plot of the image, showing different sky and extraction regions
##########

def plot_region(file, object, date, section_sky, section_spect):

# Load the FITS file
 
        hdulist = fits.open(file)
        image = hdulist['SCI'].data

# Define the size of the plot to make, in pixels

        year_obs = get_year_obs(file)

	if (year_obs == 2014):
          dy = 300 # 300 -> 150 pixels above and below the spectrum itself

	if (year_obs == 2015):
          dy = 200 # 300 -> 150 pixels above and below the spectrum itself

        dx = 300

# Define the size of the plot, in cm

        rcParams['figure.figsize'] = 10, 10 # Make plot normal size

# Calc the LHS and RHS of the image
        
        x0 = 5
        x1 = image.shape[1] - x0

# Calc the top/bottom of the regions to extract for spectrum and sky

# Extract two halves of the image, and paste them together
        
        image_extract_1 = image[:, 0:dx/2] # Left portion of image
        image_extract_2 = image[:, -dx/2:] # Right portion of image
        image_extract = np.hstack((image_extract_1, image_extract_2))

        fig = plt.imshow(log(image_extract), origin='lower', interpolation='none') # 'lower' means that row order matchs that in DS9

# Define the size of the plot, in cm

        rcParams['figure.figsize'] = 10, 10 # Make plot normal size

# Calc the LHS and RHS of the image
        
        x0 = 5
        x1 = image.shape[1] - x0

# Calc the top/bottom of the regions to extract for spectrum and sky

# Extract two halves of the image, and paste them together
        
        image_extract_1 = image[:, 0:dx/2] # Left portion of image
        image_extract_2 = image[:, -dx/2:] # Right portion of image
        image_extract = np.hstack((image_extract_1, image_extract_2))

        fig = plt.imshow(log(image_extract), origin='lower', interpolation='none') # 'lower' means that row order matchs that in DS9
         
        x1 = dx - x0

# Parse the section_sky variable and extract its values
         
        y_sky_plot   = np.array(string.split(string.replace(string.replace(section_sky,   '[', ''), ']', ''), ':'), dtype='float')
        y_spect_plot = np.array(string.split(string.replace(string.replace(section_spect, '[', ''), ']', ''), ':'), dtype='float')
        
        plot([x0, x1, x1, x0, x0], [y_spect_plot[0], y_spect_plot[0], y_spect_plot[1], y_spect_plot[1], y_spect_plot[0]], 
             color='black')
        plot([x0, x1, x1, x0, x0], [y_sky_plot[0],   y_sky_plot[0],   y_sky_plot[1],  y_sky_plot[1],    y_sky_plot[0]],   
             color='blue')
        plt.ylim((y_spect_plot[0]-dy/2, y_spect_plot[0] + dy/2))
        plt.xlim((0, dx))
   	
        plt.title(shorten_filename(file) + ", " + object + " " + date + 
 	                                   ", SKY = " + section_sky + ", SPECTRUM = " + section_spect, fontsize=15)
        plt.text(dx/2, y_sky_plot[0],   "SKY: "   + section_sky, color='white', fontsize=15)
        plt.text(dx/2, y_spect_plot[0], "SPECT: " + section_spect, color='white', fontsize=15)

        plt.show()
        
        hdulist.close()

##########
# function call_specsky
##########

def call_specsky(file, i, section_sky, file_out, verbose=0):

    if (verbose > 0):
        print "Image " + repr(i) + '/' + repr(size(files)) + ': ' + file + \
         ', section = ' + section_sky
    
        print "\n ** Calling SPECSKY(" + file + "), section = " + section_sky + "\n\n"
        
    s = iraf.specsky(images=file, outimages=file_out, section=section_sky, clobber=1,
                    logfile='salt_specsky.log', outpref='', verbose=(verbose > 1))


##########
# function call_specrectify
##########

def call_specrectify(file, i, object, calfile, file_out, verbose=0):

      if (verbose > 0):
        print "\n ** Calling SPECRECTIFY(" + file + ")\n\n"

      s = iraf.specrectify(images=file, solfile=calfile, outimages=file_out, caltype='line', 
          function='polynomial', order=3, inttype='interp', outpref='', w1='None', w2='None', dw='None', nw='None', 
          blank=0.0, clobber=1, logfile='salt_specrectify.log', verbose=(verbose > 1))

      rcParams['figure.figsize'] = 20, 30 
      hdulist = fits.open(file_out)
      image = hdulist['SCI'].data
      dy_spect_plot = 400		# How many rows vertically is the spectrum? ** Changes 2014 vs. 2015
      ny = image.shape[0]	# Look up the central row of the spectrum, and take +- dy/2 rows from that
      plt.imshow(log(image))
      plt.ylim((ny/2 - dy_spect_plot/2, ny/2 + dy_spect_plot/2))

      plt.title(file + ", " + object + ", " + ' RECT')
      plt.show()
      hdulist.close()


##########
# Funtion CALL_CRCLEAN
##########

def call_crclean(file, i, object, file_cr, mode_cr, iter_cr, verbose=0):
    
    if (object not in ["Pluto", "HD 146233", "HD146233", "ARC"]):
      print "Unknown object -- can't CRCLEAN it"
      return
    else:
      if (verbose > 0):
        print "\n ** Calling CRCLEAN(" + files[i] + "), i = " + repr(i) + ", object = " + object + "\n\n"

      s = iraf.saltcrclean(files[i], files_cr[i], '', crtype=mode_cr, maxiter = iter_cr, 
                                 logfile = 'salt_crclean.log', verbose=(verbose > 1), clobber=1)

##########
# FUNCTION CALL_SPECIDENTIFY
##########

def call_specidentify(file, i, calfile, object, grating, tilt, lampid, year_obs, dir_linelists, verbose=0):

    if (object != "ARC"):
      print "Can't do SPECIDENTIFY: Object = " + object + ", not ARC"
      return

    if (verbose > 0):
      print "\n ** Calling SPECIDENTIFY, " + file + ", i = " + repr(i) + ", object = " + object + \
      ", lampid = " + \
      lampid.replace(" ", "") + ", grating = " + grating + ", tilt = " + repr(tilt) + ' deg' + "\n\n"
 
    a = tilt
    if (tilt < 15):
      guessfile = 'salt_pluto_cal_guess_angle_13.txt'
    else:
      guessfile = 'salt_pluto_cal_guess_angle_20.txt'

    mdiff = 10 # differnce between observed line position and catalog line position. In tutorial, this was 5, 
	       # and I missed a lot of lines as a result.

# Figure out which ARC light file to use. As per Steve Crawford, I should use the .salt ones, not the .txt ones.

    if (year_obs == 2015):
      rmiddle = 342			# 2015 data seemed to be binned more aggressively than 2014, with fewer rows in the output!
      rstep = 50
      rstart = rmiddle - rstep

    out = iraf.specidentify(images=file, linelist=dir_linelists + '/' + lampid.replace(" ", "") + '.salt', 
	outfile=calfile, automethod='MatchZero',
   	guesstype='rss', guessfile=None,
#   	guesstype='file', guessfile=guessfile,
	smooth=3, 
   	function='polynomial', order=3, rstep=rstep, rstart=rstart,
#   	function='polynomial', order=3, rstep=rstep, rstart='middlerow',		# rstart='middlerow' for v0.47. In nightly can use #.
	mdiff=mdiff, thresh=3, startext=0, niter=5, inter=1, clobber=0, logfile='salt.log', verbose=(verbose > 1))

##########
# function call_specextract
##########

def call_specextract(file, i, section_spect, file_out, verbose=0):

    # Now that we have calculated the position of the spectrum, put it into the output filename
    # Spectral position is used for both sky and spectrum extraction, but we only list it here.

        if (verbose > 0):
          print "\n ** Calling SPECEXTRACT(" + file + ") , section = " + section_spect + "\n\n"      

        iraf.specextract(images=file, outfile=file_out, method='normal', 
          section=section_spect,
          thresh=3.0, minsize=3.0, 
          outformat='ascii', convert=1, clobber=True, logfile='salt_spectextract.log', verbose=(verbose > 1))

##########
# Function to detect if specidentify has been run on a file
##########

def exists_specidentify(file, calfile, object):
  if (object != 'ARC'):
    return -1
  return  (subprocess.call(['grep', file, calfile]) == 0)	# Has file been run thru SPECIDENTIFY yet?

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
        
dir_data_2014 = "/Users/throop/Data/SALT_Pluto_2014/product"
dir_data_2015 = "/Users/throop/Data/SALT_Pluto_2015/product"

dir           = '/Users/throop/python'
dir_pysalt    = '/Users/throop/iraf/pysalt'
dir_linelists = dir_pysalt + "/data/linelists"

##########
# Set flags for which of the tasks we want to do when we hit 'Do All'. These flags 
# can be set and cleared from within the program
##########

DO_SPECIDENTIFY = False
DO_CRCLEAN      = True
DO_SPECRECTIFY  = True
DO_SPECSKY      = True
DO_SPECEXTRACT  = True
DO_PLOT_IMAGE   = True # Make a plot showing the output image, with sections marked, etc
DO_PLOT_RECTIFY = True # Make a plot showing the rectified spectrum alone, so we can see how good the SPECIDENTIFY was
DO_PLOT_SPECTRUM= True # Plot the output spectrum itself

DO_ARCS_ONLY    = True  # Downselect to only ARCs
DO_ALL_BUT_ARCS = False   # Downselect to only science frames

verbose		= 1	# 0: no messages
			# 1: Messages from salt_calibrate, but not PySALT
			# 2: Messages from salt_calibrate and PySALT
##########
# Get some input from the user
##########

inp = raw_input("Year: 201(4) or 201(5): ")
if (inp.find('4') >= 0): year_obs = 2014
if (inp.find('5') >= 0): year_obs = 2015

inp = raw_input("(a)rcs or (d)ata: ")
if (inp == 'a'):
  DO_ARCS_ONLY    = True
  DO_ALL_BUT_ARCS = False

if (inp == 'd'):
  DO_ARCS_ONLY    = False
  DO_ALL_BUT_ARCS = True

##########
# Set Parameters for all of the steps
##########

#   year_obs        = 2015  # Set this to 2014 or 2015, depending on which dataset

mode_cr         = 'median'		# Also can be 'edge' or 'fast', but I think that median is best
iter_cr         = 3
y0_sky_default  = 50      # Position for sky spectral extraction start. 0 = center of Pluto / HD spectrum. Can be overridden.

dy_sky          = 20      # Number of rows to use for the sky spectrum
dy_spect        = 20      # Number of rows to use for the spectrum itself (full width)

if (year_obs == 2015):    # Set the offset position to start the spectrum. I think different for 2014 vs 2015 because of binning?
  y0_sky_default = 25

if (year_obs == 2014):
  y0_sky_default = 50


# CD into python. Just dump everything there for now.
# I don't know how PyRAF handles directories, so better to cd into the right one from the start.

os.chdir(dir)

# Create a list for each of the keywords. After we do the loop, 
# each list will be populated -- e.g., FITS_EXPTIME = [1, 1, 10, 1], etc

fits_exptime = []	# new list (same idea as array) 
fits_object  = [] 
fits_date_obs= []
fits_utc_obs = [] # starting time of exposure
fits_jd      = [] 
fits_obsmode = [] 
fits_moonang = []  
fits_instrume= [] 
fits_gainset = [] 
fits_rospeed = [] 
fits_ccdsum  = [] # '2 2' -> 2x2 bininng 
fits_pixscale= [] # before binning 
fits_proposer= [] 
fits_propid  = [] 
fits_airmass = [] 
fits_masktyp = [] 
fits_maskid  = []
 
fits_ra      = [] 
fits_dec     = [] 
fits_telra   = [] 
fits_teldec  = [] 

fits_grating = [] 
fits_gr_sta  = [] 
fits_gr_angle= []
fits_lampid=   []

# Get a list of raw FITS files. Do not include *_rect_sky.fits, etc -- just the 'original' fits files.
# Also, it's fine to exclude FLATs, since we can basically ignore them entirely in the pipeline as per Steve Crawford 26-Aug-14.

if (year_obs == 2014):
  file_list = glob.glob(dir_data_2014 + '/mbxgpP2014[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits') # glob won't take regexp

if (year_obs == 2015):
  file_list = glob.glob(dir_data_2015 + '/mbxgpP2015[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits') # glob won't take regexp

files = np.array(file_list, dtype='S200') # Need to make string fields long enough for expansion!

# Read the JD from each file. Then sort the files based on JD.

# files      = np.array(file_list,dtype='S100')
 
print "Reading " + repr(size(file_list)) + " files"
jd = []
object = []
for file in files:
    sys.stdout.write('.') # print w/o newline
    hdulist = fits.open(file)
    jd.append(hdulist[0].header['JD']) # For now, only read JD and OBJECT
    object.append(hdulist[0].header['OBJECT'])

    hdulist.close()
    
print
print "Sorting by JD, removing flats, removing bad files"

object = np.array(object, dtype='S30')

#   quit

if (DO_ARCS_ONLY):
    arcs = [object == 'ARC']
    files = files[arcs]
    object = object[arcs]
    jd = np.array(jd)
    jd = jd[arcs]

    print
    print 'Doing ARCs only...'

if (DO_ALL_BUT_ARCS):
    arcs = [object != 'ARC']
    files = files[arcs]
    object = object[arcs]
    jd = np.array(jd)
    jd = jd[arcs]

    print
    print 'Doing all but ARCs...'

indices = np.argsort(jd)
files = files[indices]          # Resort by JD    
files = files[object != 'FLAT'].copy() # Remove all of the flats

# Remove some bad files on a one-off basis.
# NB: My philosophy here is to remove as little as possible. Process what I can. Then in the next routine, I can 
#     mix & match.

# files = files[files != dir_data + '/mbxgpP201409050032.fits'] # Pluto spectra. Blank. 0029, 0030, 0031 are all OK.

# 20140905 0{189, 190, 191} : All three of the these are on-target but very faint.
# 20140905 0{198, 199, 200, 201}: All four of these are usable. Not sure why we have four. 198, 199 = best of the bunch.

# files = files[files != dir_data + '/mbxgpP201409050189.fits'] # Pluto spectra. Took four. First two are clearly off-target
# files = files[files != dir_data + '/mbxgpP201409050190.fits'] # Pluto spectra. Took four. First two are clearly off-target
# files = files[files != dir_data + '/mbxgpP201409050198.fits'] # Very faint. Off target? 200 / 201 are OK. Took four.
# files = files[files != dir_data + '/mbxgpP201409050199.fits'] # Very faint. Off target? 200 / 201 are OK. Took four.

# *** 191 looks off-target also??

print

print "Reading " + repr(size(files)) + " files"

for file in files:
    sys.stdout.write('.')
        
    hdulist = fits.open(file)
    
    header = hdulist[0].header
    keys = header.keys()

    fits_object.append(header['OBJECT'])
    fits_exptime.append(header['EXPTIME'])
    fits_date_obs.append(header['DATE-OBS'])
    fits_utc_obs.append(header['UTC-OBS'])
    fits_jd.append(header['JD'])
    fits_obsmode.append(header['OBSMODE'])
    fits_moonang.append(header['MOONANG'])
    fits_instrume.append(header['INSTRUME'])
    fits_gainset.append(header['GAINSET'])
    fits_rospeed.append(header['ROSPEED'])
    fits_ccdsum.append(header['CCDSUM'])
    fits_pixscale.append(header['PIXSCALE'])
    fits_proposer.append(header['PROPOSER'])
    fits_propid.append(header['PROPID'])
    fits_airmass.append(header['AIRMASS'])
    
    fits_dec.append(header['RA'])
    fits_ra.append(header['DEC'])
    fits_telra.append(header['TELRA'])
    fits_teldec.append(header['TELDEC'])
    
    if ('GRATING' in keys):    
        fits_grating.append(header['GRATING'])
        fits_gr_sta.append(header['GR-STA'])
        fits_gr_angle.append(header['GR-ANGLE'])
        fits_masktyp.append(header['MASKTYP']) # not always there
        fits_maskid.append(header['MASKID'])
    else:
       fits_grating.append('')
       fits_gr_sta.append('')
       fits_gr_angle.append('')
       fits_masktyp.append('')
       fits_maskid.append('')

    if ('LAMPID' in keys):
        fits_lampid.append(header['LAMPID'])
    else:
        fits_lampid.append('')

# Convert from lists to NumPy arrays as needed

jd          = np.array(fits_jd)
utc         = np.array(fits_utc_obs)
date        = np.array(fits_date_obs)
airmass     = np.array(fits_airmass)
object      = np.array(fits_object, dtype='S30') # np.array can use string arrays as easily as float arrays
instrument  = np.array(fits_instrume)
proposer    = np.array(fits_proposer)
exptime     = np.array(fits_exptime)
lampid	  = np.array(fits_lampid)
grating     = np.array(fits_grating)
tilt        = np.array(fits_gr_angle)
obsmode     = np.array(fits_obsmode)  # Remove extra spaces

# Extract just the good ones
# ** When initializing the np string array, need to set the maximum string length.
# ** I choose a long default here, so we won't exceed it.

# files      = np.array(file_list,dtype='S100')

files_short = np.array(files)

for i in range(files.size):
    files_short[i] = files[i].split('/')[-1]  # Get just the filename itself
   
is_bad = ((proposer != 'Throop') | 
          (instrument != 'RSS') | 
          (object == 'BIAS') | 
          (obsmode == 'FABRY-PEROT') | 
          (files_short == "mbxgpP201408260016.fits"))  # This file is an ARC but it is bad... all noise. 
                                                       # Don't know why. It was retaken, so skip it.
          
is_good = (is_bad == False)

# pdb.set_trace()

files      = files[is_good]  # Is this legal in python? Does it make sense, or is it a view into itself??
files_short=files_short[is_good]
jd         = jd[is_good]
utc        = utc[is_good]
date       = date[is_good]
airmass    = airmass[is_good]
object     = object[is_good]
instrument = instrument[is_good]
proposer   = proposer[is_good]
exptime    = exptime[is_good]
lampid     = lampid[is_good]
grating    = grating[is_good]
obsmode    = obsmode[is_good]
tilt       = tilt[is_good].astype(float) # For some reason this reads as a string from the FITS file...  

# 2014 uses 'HD 146233', while 2015 is 'HD146233'. Probably my error in making Step-2? Make them match.

for i in range(size(object)):
  oi = object[i]
  if oi == "HD146233": object[i] = "HD 146233"

# Make arrays for the various output files

print "\ndone reading files\n"

# Get a list of Arcs, Flats, etc. Every image is exactly one of these.

is_arc   = (object == 'ARC')
is_pluto = (object == 'Pluto')
is_hd    = (object == 'HD 146233')
is_flat  = (object == 'FLAT')
                 
## Call specidentify to identify individual lines in the cal spectrum

# "SPECIDENTIFY identies Arc lines in a spectrum and calculates the wavlength
# identification of the lines or it can be run automatically where the task will
# identify the lines without any assistance from the user."
#
# Q: Where do I get the line lists from A:  SAAO-SALT webpage. 

# Display a file
# iraf.imexamine(files[2] + '[1]')  # display an arc spectrum in a ds9 window.
# for file in files do:
# iraf.imhead(filesarr[is_arc][2]+'[0]',long=1)

calfile         = "salt_pluto_cal.txt" # File which gets written by SPECIDENTIFY with the parameters of the wavelength solution.
                                       # This is a running file, which gets appended every time SPECRECTIFY is run. It has many entries.
file_positions  = "salt_positions_extract.csv" # File which I made from Excel and ds9, of which rows to extract.

# Read it into an array of lists
# Problem is that numpy arrays have a single datatype

# f = loadtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, 
#            dtype = {'names' : ('djd', 'utc', 'object', 'exptime', 'tilt', 'file', 'sky', 'spectrum', 'notes'),
#                     'formats': ('f4', 'S10', 'S10',     'f4',     'f4',    'S10', 'S10', 'S10',      'S10')} )

# Read it in as a string array. Then extract the proper columns, and use those. Not very slick.

# f = loadtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, 
#            dtype = 'S300')
# genfromtxt() is more general than loadtxt()

f = genfromtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, dtype = 'str', invalid_raise=False, usecols=range(9))

name_files_positions  = f[:,5]
positions_extract_sky     = f[:,6] # Which rows of each image do I use for the sky subtraction?
positions_extract_spect   = f[:,7] # Which rows of each image do I use for the spectrum extraction?

#########
# Create the array for the sky y-positions 
#########

y0_sky          = y0_sky_default + np.zeros(size(files))

# Change these one-off values for the Y position of the sky. This is because there are stars in these positions.

# 2014 values

y0_sky[ np.char.find(files, '06280032') > 0 ] = 40
y0_sky[ np.char.find(files, '06280033') > 0 ] = 40
y0_sky[ np.char.find(files, '06280034') > 0 ] = 40
y0_sky[ np.char.find(files, '06280041') > 0 ] = 40
y0_sky[ np.char.find(files, '06280042') > 0 ] = 40
y0_sky[ np.char.find(files, '06280054') > 0 ] = 30
y0_sky[ np.char.find(files, '06280055') > 0 ] = 30

y0_sky[ np.char.find(files, '07020053') > 0 ] = 25
y0_sky[ np.char.find(files, '07020054') > 0 ] = 25
y0_sky[ np.char.find(files, '07020061') > 0 ] = 25
y0_sky[ np.char.find(files, '07020062') > 0 ] = 25
y0_sky[ np.char.find(files, '09050191') > 0 ] = 40

# 2015 values

y0_sky[ np.char.find(files, '05030054') > 0 ] = 35
y0_sky[ np.char.find(files, '05030053') > 0 ] = 35
y0_sky[ np.char.find(files, '05030051') > 0 ] = 35
y0_sky[ np.char.find(files, '05030050') > 0 ] = 35

##########
# Create filenames for all the other files we might load
##########

files_cr              = np.empty_like(files, dtype = 'S200')
files_cr_rect         = np.empty_like(files, dtype = 'S200')
files_cr_rect_sky     = np.empty_like(files, dtype = 'S200')
files_cr_rect_sky_ext = np.empty_like(files, dtype = 'S200')

# Create strings like '[333:355]' for the region of sky to extract

section_sky           = np.empty_like(files)
section_spect         = np.empty_like(files)

for i in range(files.size):

# Create the CR (cosmic-ray) filename

    files_cr[i]              = string.replace(files[i], '.fits', '_crmode=' + mode_cr + '_criter=' + repr(iter_cr) + '.fits')
    files_cr_rect[i]         = string.replace(files_cr[i], '.fits', '_rect.fits')
    files_cr_rect_sky[i]     = string.replace(files_cr_rect[i], '.fits', '_y0sky=' + repr(int(y0_sky[i])) + '_dysky=' + repr(dy_sky) + 
                                               '.fits') 

# Create the _ext.txt filename
# Q: What does this get_row_max stuff try to do??
# A: You give it a filename, and it searches the file for the brightest row, which it uses to extract the spectrum from.

    f = string.replace(files_cr_rect_sky[i], '.fits', '_dyspect=' + repr(dy_spect) + 
                                               '_y0spect=XYZ' + '.txt')
    files_cr_rect_sky_ext[i] = string.replace(f, 'XYZ', repr(get_row_max(files_cr_rect_sky[i])))

    section_sky[i]           = ''
    section_spect[i]         = ''
   
# Order of operation:
#  1. Spec identify (interactive). No parameters.
#  2. CR rejection / cleaning. Parameters: crmode=edge | fast | median; maxiter=3
#  3. Spec rectify. No parameters. (Based on fit in spec identify)
#  4. Sky subtraction
#  5. Spectral extraction
#  6. Make plot
#
# Example filename:
# image_cr.edge.iter5_rect_sky.50.10_extract

##########    
# Now start the main loop
##########

i = 0
tableheader = "#    File                   Type     Exptime SPCIDNTFY _cr  _rect  _sky  _ext   #"

print tableheader

IS_DONE = False

range_selected = range(size(files))			# By default, set the range to act upon every file

prompt = "(#)number, \n" + \
      'spec(i)dentify, (r)CRCLEAN, (t)RECTIFY, specs(k)y, spece(x)tract, \n' + \
      '(pr)plot rectify, (pi)plot image, (ps)plot spectrum, ' + \
      'perform (A)ll flagged actions, (sr)set range of files for action, (f)lags, \n' + \
      "(q)uit, (l)ist, (n)ext, (p)revious, (v)erbose,  (" + repr(i) + ") ?: "

for j in range(size(files)):
      print_status_line_file(
	      j, files_short[j], object[j], int(exptime[j]), 
	      exists_specidentify(files_short[j], calfile, object[j]),	# Has file been run thru SPECIDENTIFY?
	      os.path.isfile(files_cr[j]), 					# Does _cr.fits file exist?
	      os.path.isfile(files_cr_rect[j]),				# Does _cr_rect.fits file exist?
	      os.path.isfile(files_cr_rect_sky[j]),				# Does _cr_rect_sky.fits file exist?
	      os.path.isfile(files_cr_rect_sky_ext[j])) 			# Does _cr_rect_sky_ext.txt file exist?

while (IS_DONE == False): 

	print tableheader
	print_status_line_file(
	          i, files_short[i], object[i], int(exptime[i]), 
		  exists_specidentify(files_short[i], calfile, object[i]),
		  os.path.isfile(files_cr[i]), 
		  os.path.isfile(files_cr_rect[i]),
		  os.path.isfile(files_cr_rect_sky[i]),
		  os.path.isfile(files_cr_rect_sky_ext[i]))
       
        print

	inp = raw_input("File " + repr(i) + ":" + prompt)

        if (inp == 'sr'):
       
	    inp2 = raw_input("Range of files [e.g. *; 1-10; 44]: ")
	    if (inp2.find('-') >= 0):		# Range of files
		min = int((inp2.split('-'))[0])
		max = int((inp2.split('-'))[1])
		range_selected = range(min, max+1)

	    elif (inp2.find('*') >= 0):		# All files
		range_selected = range(size(files))

	    else:		# Only a single file
		range_selected = [int(inp2)]

	    print 'Range: ' + repr(range_selected)
		
	if (inp == 'r'):
            call_crclean(files[i], i, object[i], files_cr[i], mode_cr, iter_cr, verbose=verbose)
	  
	if (inp == 'i'):
            call_specidentify(files[i], i, calfile, object[i], grating[i], tilt[i], lampid[i], year_obs, dir_linelists, verbose=verbose)

	if (inp == 't'):
            call_specrectify(files_cr[i], i, object[i], calfile, files_cr_rect[i], verbose=verbose)
	  
	if (inp == 'k'):
	    section_sky = get_section_sky(files_cr_rect[i], y0_sky[i], dy_sky) # returns a string, like "[1:100]"
            call_specsky(files_cr_rect[i], i, section_sky, files_cr_rect_sky[i], verbose=verbose)

	if (inp == 'x'):
#               files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky_ext[i], 'XXX', repr(get_row_max(files_cr_rect_sky[i])))
	    section_spect = get_section_spect(files_cr_rect[i], dy_spect)
	    print "About to extract using section " + repr(section_spect)
            call_specextract(files_cr_rect_sky[i], i, section_spect, files_cr_rect_sky_ext[i], verbose=verbose)
	  
        if (inp == 'q'):
            quit()
        
	if (inp == 'f'):    # Set flags
	  print DO_SPECIDENTIFY*1, DO_CRCLEAN*1, DO_SPECRECTIFY*1, DO_SPECSKY*1, DO_SPECEXTRACT*1, \
	        DO_PLOT_RECTIFY*1, DO_PLOT_IMAGE*1, DO_PLOT_SPECTRUM*1
          print "SPECIDENTIFY, CRCLEAN, SPECRECTIFY, SPECSKY, SPECEXTRACT, PLOT_RECTIFY, PLOT_IMAGE, PLOT_SPECTRUM"
          inp2 = raw_input("Enter flag pattern: ")
	  DO_SPECIDENTIFY  = (inp2[0] == '1')
	  DO_CRCLEAN       = (inp2[1] == '1')
	  DO_SPECRECTIFY   = (inp2[2] == '1')
	  DO_SPECSKY       = (inp2[3] == '1')
	  DO_SPECEXTRACT   = (inp2[4] == '1')
	  DO_PLOT_RECTIFY  = (inp2[5] == '1')
	  DO_PLOT_IMAGE    = (inp2[6] == '1')
	  DO_PLOT_SPECTRUM = (inp2[7] == '1')
	  print DO_SPECIDENTIFY*1, DO_CRCLEAN*1, DO_SPECRECTIFY*1, DO_SPECSKY*1, DO_SPECEXTRACT*1, \
	        DO_PLOT_RECTIFY*1, DO_PLOT_IMAGE*1, DO_PLOT_SPECTRUM*1

        if (inp == 'v'):
	  inp = raw_input("Verbosity level: ")
	  verbose = int(inp)
	  print "Verbosity changed to " + repr(inp)

        if (inp == 'pi'):				# Plot image
	    section_sky   = get_section_sky(files_cr_rect[i], y0_sky[i], dy_sky)
	    section_spect = get_section_spect(files_cr_rect[i], dy_spect)
            plot_region(files_cr_rect[i], object[i], date[i], section_sky, section_spect)

        if (inp == 'ps'):				# Plot spectrum
#               files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky_ext[i], 'XXX', repr(get_row_max(files_cr_rect_sky[i])))
            plot_spectrum(files_cr_rect_sky_ext[i], object[i], date[i])

        if (inp == 'pr'):				# Plot spectrum
#               files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky_ext[i], 'XXX', repr(get_row_max(files_cr_rect_sky[i])))
            plot_rectify(files_cr_rect[i], object[i], date[i])

        if (inp == 'n'):			# Next file
	  i+=1

        if (inp == 'p'):			# Previous file
	  i-=1

	if (inp == 'A'):		# Process all files, according to flags set
	  prompt2 = ' Applying all flagged steps to ' + repr(size(range_selected)) + ' selected files. Proceed? (n)'
	  inp2 = raw_input(prompt2)
	  if (inp2 == 'y'):
	    for j in range_selected:
#   	    range(size(files)):

		if (DO_CRCLEAN):      
		  call_crclean(files[j], j, object[j], files_cr[j], mode_cr, iter_cr, verbose=verbose)

		if (DO_SPECIDENTIFY): 
		  call_specidentify(files[j], j, calfile, object[j], grating[j], tilt[j], lampid[j], year_obs, dir_linelists, verbose=verbose)

		if (DO_SPECRECTIFY):  
		  call_specrectify(files_cr[j], j, object[j], calfile, files_cr_rect[j], verbose=verbose)

		if (DO_SPECSKY): 
		  section_sky = get_section_sky(files_cr_rect[j], y0_sky[j], dy_sky) # returns a string, like "[1:100]"
		  call_specsky(files_cr_rect[j], j, section_sky, files_cr_rect_sky[j], verbose=verbose)

		if (DO_SPECEXTRACT):
#                     files_cr_rect_sky_ext[j] = string.replace(files_cr_rect_sky_ext[j], 'XXX', repr(get_row_max(files_cr_rect_sky[j])))
	          section_spect = get_section_spect(files_cr_rect[j], dy_spect)
		  call_specextract(files_cr_rect_sky[j], j, section_spect, files_cr_rect_sky_ext[j], verbose=verbose)

                if (DO_PLOT_IMAGE):
	          section_sky   = get_section_sky(files_cr_rect[j], y0_sky[j], dy_sky)
	          section_spect = get_section_spect(files_cr_rect[j], dy_spect)
                  plot_region(files_cr_rect[j], object[j], date[j], section_sky, section_spect)

                if (DO_PLOT_RECTIFY):
		  plot_rectify(files_cr_rect[j], object[j], date[j])
		  
   		if (DO_PLOT_SPECTRUM):
                  plot_spectrum(files_cr_rect_sky_ext[j], object[j], date[j])

        if ((inp == 'l') | (inp == '?')):
            print tableheader
            for j in range(size(files)):

              print_status_line_file(
	          j, files_short[j], object[j], int(exptime[j]), 
		  exists_specidentify(files_short[j], calfile, object[j]),
		  os.path.isfile(files_cr[j]), 				
		  os.path.isfile(files_cr_rect[j]),		
		  os.path.isfile(files_cr_rect_sky[j]),	
		  os.path.isfile(files_cr_rect_sky_ext[j])) 

        if (is_number(inp)):
	    i = int(inp)
	    print "n = " + repr(i)

##########
# Step 1: SPECIDENTIFY. Do this for ARCs *only*.
#
# This is the one where we identify spectral lines, and do the actual wavelength calibration
#
# How to use SPECIDENTIFY: Hit 'Arc', then 'a', 'b', 'f', 'X'. Residual. 'q'. Residuals are typically < 1.
##########

#       if (DO_SPECIDENTIFY):
#         if (exists_specidentify(files[i], calfile)):
#           print "Already exists! Process again?"
#         else:
#           call_specidentify(files[i], i, calfile, object[i], grating[i], tilt[i], lampid[i], year_obs, dir_linelists)

##########
# Step 2: Process cosmic rays. 
# Not really necessary on ARCs, but we do it just for consistency to keep filenames the same.
##########

#       if (DO_CRCLEAN):
#         call_crclean(files[i], i, object[i], files_cr[i], iter_cr)

##########
# Step 3: Spec Rectify. Call this for all SCIENCE data (HD and Pluto).
# Can also do on ARCs. It is not necessary, but it shows us how well the SPECIDENTIFY worked.
# This uses the polynomial transformations above to rewrite new FITS files which are wavelength calibrated.
# Seems like a backwards way to do it (aren't there artifacts introduced here?), but that's the process, so OK.
# ** Steve's walkthru used function=polynomial above and function=legendre below. That is an error.
#
# This routine *automatically* finds the proper solution, from the 
#    "The task will find the calibration using the same instrument setup that is closest in time 
#     and use that calibration to correct the data."
#
# This just runs. No input required. Takes about 30 sec per image frame.
# I did get an error on one file: "Improper input: N=4 must not exceed M=0"
# ** Problem image: mbxgpP201406240026.fits: "TypeError: Improper input: N=4 must not exceed M=0" in minpack.pyc
# ** Solved: problem was that this was a Fabrey-Perot image. That was an error to take them.
##########

#       if (DO_SPECRECTIFY):
#         call_specrectify(files[i], i, object[i], calfile, files_cr_rect[i])

##########       
# Step 4: SPECSKY: Remove the sky from a portion of the image.
# Apply this to science data (ARC and Pluto). 
# Looks like this *cannot* take an input file which is a list of 'section's to use for the sky subtraction.
# I need to manage that myself.
# *** To do: Read the 'section' values from a file. They should be computed for every image.
##########

#       if (DO_SPECSKY):
#         call_specsky(files[i], i, file_out)
    
##########
# SPECEXTRACT: Now that everything is finished, extract the spectrum and write to txt file
# *** To do: read the 'section' values from a file. They should be computed for every image.
##########

#       if (DO_SPECEXTRACT): 
#         row_max = get_row_max(files_cr_rect_sky_ext[i])
#         files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky_ext[i], 'XXX', row_max)
#         section_spect = get_section_spect(file, dy_spect)
#         call_specexctract(files_cr_rect_sky[i], i, section_spect, files_cr_rect_sky_ext[i])

##########
# Step 6: Make a plot of the image, showing where the sky and spectral extraction regions are
##########

#       section_sky   = get_section_sky(files_cr_rect[i], y0_sky[i], dy_sky)
#       section_spect = get_section_spect(files_cr_rect[i], dy_spect)
#       plot_region

##########
# Now plot the spectra
##########

files_spec = files_cr_rect_sky_ext # We've already created the filenames

fluxes        = []
dfluxes       = []
wavelengths   = []

for i in range(size(files_spec)):
    file = files_spec[i]
    if (object[i] == 'ARC'):
        print "Skipping ARC: " + file
        wavelengths.append([])
        fluxes.append([])
        dfluxes.append([])
        
    else: 
        print "Reading spectrum " + file
        d = loadtxt(file)
        wavelengths.append(d[:,0])
        fluxes.append(d[:,1])
        dfluxes.append(d[:,2])

rcParams['figure.figsize'] = 20, 30 # Make plot normal size

colors = np.array(['grey', 'brown', 'orange', 'yellow', 'blue', 'red', 'pink', 'green'])
for i in range(len(fluxes)):
    if (object[i] == 'Pluto'):
        f = np.array(fluxes[i]).astype(float)
        f = f / np.amax(f)  # normalize it
        plot(wavelengths[i],f+i, color=colors[i % size(colors)], ls='-')
    plt.axis([3000,9000,0,100])
  
plt.xlim((3000,9000))
        
##########
# Now extract the HD spectra
##########
        
w_spect_hd_blue  = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 13) & (tilt < 14))[0]
w_spect_hd_red   = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 20) & (tilt < 21))[0]

w_spect_pluto_blue  = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 13) & (tilt < 14))[0]
w_spect_pluto_red   = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 20) & (tilt < 21))[0]

