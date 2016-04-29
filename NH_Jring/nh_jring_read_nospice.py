# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

@author: throop
"""

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
import math       # We use this to get pi. Documentation says math is 'always available' but apparently it still must be imported.
from   subprocess import call
import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
#import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import wcsaxes
import HBT as hbt
import warnings
import pdb
import imreg_dft as ird

#function calc_median_framelist(t)

def get_jring_points_radec(et, num_pts=100, radius = 122000):
    "Get an array of points of RA, Dec for the Jupiter ring, at a specified radius, seen from NH at the given ET."
    
# Now calculate the ring points...

    radii_ring = np.array([1220000., 129000.])  # Inner and outer radius to plot
    num_pts_ring = 100
    
    ring_lon = np.linspace(0, 2. * np.pi, num_pts_ring)
    ra_ring  = np.zeros(num_pts_ring)
    dec_ring = np.zeros(num_pts_ring)
    
    frame = 'J2000'
    abcorr = 'LT'
#    rot = cspice.pxform('IAU_Jupiter', frame, et) # Get matrix from arg1 to arg2
    rot = [[1,0,0],[0,1,0],[0,0,1]]
    
#    st,ltime = cspice.spkezr('Jupiter', et, frame, abcorr, 'New Horizons')
    st = [0,0,0,0,0,0]
    ltime = 0
    pos = st[0:3]
    vel = st[3:6] # velocity, km/sec, of jupiter
    
    for radius_ring in radii_ring:
        for j in range(num_pts_ring):
            xyz = np.zeros(3)
            xyz = np.array((radius_ring * np.cos(ring_lon[j]), radius_ring * np.sin(ring_lon[j]), 0.))
            
            d_j2000_xyz = np.dot(rot,xyz)  # Manual says that this is indeed matrix multiplication
            j2000_xyz = 0 * d_j2000_xyz
            j2000_xyz[0:3] = d_j2000_xyz[0:3] # Looks like this is just copying it

            rho_planet    = pos                     # Position of planet
            rho_ring      = rho_planet + j2000_xyz  # Vector obs-ring
#            dist_ring     = cspice.vnorm(rho_ring)*1000 # Convert to km... CHECK UNITS!
            dist_ring = 0
            
#            range_out, ra, dec = cspice.recrad(rho_ring) # 'range' is a protected keyword in python!
            range_out = 1
            ra = 1
            dec = 1
            ra_ring[j] = ra     # save RA, Dec as radians
            dec_ring[j] = dec
                
    return ra_ring, dec_ring
    
def get_range(maxrange = 10000):
    "Request a range of input values from the user. No error checking."
    "  *; 1-10; 44   are all valid inputs. If  *  then output is 1 .. maxrange."
    
    inp2 = raw_input("Range of files [e.g. *; 1-10; 44]: ")

    if (inp2.find('-') >= 0):        # Range of files
        min = int((inp2.split('-'))[0])
        max = int((inp2.split('-'))[1])
        range_selected = range(min, max+1)

    elif (inp2.find('*') >= 0):        # All files
        range_selected = range(maxrange)

    else:        # Only a single file
        range_selected = [int(inp2)]
        print 'Range: ' + repr(range_selected)
            
    return range_selected

def match_image_brightness(arr1, arr2, frac=0.5):
    "Performs linear regression on two images to try to match them."
    "Returns fit parameter r: for best fit, use arr2 * r[0] + r[1]"
    arr1_filter = hbt.remove_brightest(arr1, frac) # Remove the stars, rings, stray light... almost anything except bg
    arr2_filter = hbt.remove_brightest(arr2, frac)
    r = linregress(arr1_filter.flatten(), arr2_filter.flatten())
    
    m = r[0] # Multiplier = slope
    b = r[1] # Offset = intercept

    arr2_fixed = arr2 * m + b
    return (m,b)
    
def get_image_nh(file, frac=0.9, polyfit=True, raw=False, ):
    " Reads an NH FITS file from disk. Does simple image processing on it, removes background,"
    " and roughly scales is for display."
    " Might still need a log scaling applied for J-ring images."
       
#        if ('hdulist') in locals():        # If there is already an hdulist, then close it. (Might double-close it, but that is OK.)
#            hdulist.close()
#        f = t['Filename'][i]   
    hdulist = fits.open(file)
    image = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'

# If requested, return the raw unscaled image

    if (raw == True):
        return arr
               
# Clip any pixels brighter than this fractional amount (e.g., 0.9 = clip to 90th percentile)
               
    arr_filtered = hbt.remove_brightest(image.data, frac)

    if (polyfit == False):
        return arr_filtered
        
# Fit the data using a polynomial

    if hdulist[0].header['TARGET'] == 'IO':
        power = 0               # Only do the polynomial fit for ring (Target = Jupiter)
    else:
        power = 10       

    polyfit = hbt.sfit(arr_filtered, power)
        
    return arr_filtered - polyfit
        
def find_stars(im):
    "Locate stars in an array, using DAOphot. Returns N x 2 array with xy positions. No magnitudes."
         
    mean, median, std = sigma_clipped_stats(im, sigma=3.0, iters=5)
    sources = daofind(im - median, fwhm=3.0, threshold=5.*std)
    x_phot = sources['xcentroid']
    y_phot = sources['ycentroid']
        
    points_phot = np.transpose((x_phot, y_phot)) # Create an array N x 2

    return points_phot
    
def calc_offset_points(points_1, points_2, shape, plot=False):
    "Calculate the offset between a pair of ordered points -- e.g., an xy list of star positions, and and xy list of model postns."
    "Returned offset is integer pixels as tuple (dy, dx)."
    
    diam_kernel = 5 # If this is 11, that is too big, and we gt the wrong answer. Very sensitive.

    image_1 = hbt.image_from_list_points(points_1, shape, diam_kernel)
    image_2 = hbt.image_from_list_points(points_2, shape, diam_kernel)
 
    t0,t1 = ird.translation(image_1, image_2) # Return shift, with t0 = (dy, dx). t1 is a flag or quality or something.
    (dy,dx) = t0
    
    if (plot):

        xrange = (0, shape[0])
        yrange = (0, shape[1])

        figs = plt.figure()
        ax1 = figs.add_subplot(1,2,1) # nrows, ncols, plotnum. Returns an 'axis'
        ax1.set_aspect('equal') # Need to explicitly set aspect ratio here, otherwise in a multi-plot, it will be rectangular

#        fig1 = plt.imshow(np.log(image_1))
        plt.plot(points_1[:,0], points_1[:,1], marker='o', color='lightgreen', markersize=4, ls='None', label = 'Photometric')
        plt.plot(points_2[:,0], points_2[:,1], marker='o', color='red', markersize=4, ls='None', label = 'Cat')
        plt.legend()
       
        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        plt.title('Raw')
        
        ax2 = figs.add_subplot(1,2,2) # nrows, ncols, plotnum. Returns an 'axis'
        plt.plot(points_1[:,0], points_1[:,1], marker='o', color='lightgreen', markersize=9, ls='None')
        plt.plot(points_2[:,0] + t0[1], points_2[:,1] + t0[0], marker='o', color='red', markersize=4, ls='None')
        ax2.set_aspect('equal')

        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        plt.title('Shifted, dx=' + repr(dx) + ', dy = ' + repr(dy))
        
        plt.show()
        
    return t0
        
hbt.set_plot_defaults()
  
d2r = np.pi /180.
r2d = 1. / d2r
         
dir_data = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'
# Start up SPICE

file_tm  = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"
#cspice.furnsh(file_tm)

# Get the full list of files

file_list = glob.glob(dir_data + '/*fit')
files = np.array(file_list)
indices = np.argsort(file_list)
files = files[indices]

DO_FAST = True
if (DO_FAST):
    files = files[0:100]
    
# Read the JD from each file. Then sort the files based on JD.

jd = []
for file in files:
    hdulist = fits.open(file)
    jd.append(hdulist[0].header['MET'])
    hdulist.close()
     
fits_met = []    # new list (same as array) 
fits_startmet  = [] 
fits_stopmet= []
fits_exptime = [] # starting time of exposure
fits_target  = [] 
fits_reqdesc = [] 
fits_spcinst0= [] 
fits_spcutcjd= []   
fits_naxis1= [] 
fits_naxis2 = [] 
fits_spctscx = [] # sc - target, dx 
fits_spctscy = [] # dy
fits_spctscz = [] # dz
fits_spctcb  = [] # target name
fits_spctnaz = [] # Pole angle between target and instrument (i.e., boresight rotation angle)

#files_short = np.array(files)
#for i in range(files.size):
#    files_short = files[i].split('/')[-1]  # Get just the filename itself

# Set up one iteration variable so we don't need to create it over and over
num_obs = np.size(files)
i_obs = np.arange(num_obs)

for file in files:
    print "Reading file " + file

    hdulist = fits.open(file)
    header = hdulist[0].header
    
    keys = header.keys()

    fits_met.append(header['MET'])
    fits_exptime.append(header['EXPTIME'])
    fits_startmet.append(header['STARTMET'])
    fits_stopmet.append(header['STOPMET'])
    fits_target.append(header['TARGET'])
    fits_reqdesc.append(header['REQDESC'])
    fits_spcinst0.append(header['SPCINST0'])
    fits_spcutcjd.append( (header['SPCUTCJD'])[3:]) # Remove the 'JD ' from before number
    fits_naxis1.append(header['NAXIS1'])
    fits_naxis2.append(header['NAXIS2'])
    fits_spctscx.append(header['SPCTSCX'])
    fits_spctscy.append(header['SPCTSCY'])
    fits_spctscz.append(header['SPCTSCZ'])    
    fits_spctnaz.append(header['SPCTNAZ'])    
       
    hdulist.close() # Close the FITS file

print object
print "done"

# Calculate distance to Jupiter in each of these
# Calc phase angle (to Jupiter)
# Eventually build backplanes: phase, RA/Dec, etc.
# Eventually Superimpose a ring on top of these
#  ** Not too hard. I already have a routine to create RA/Dec of ring borders.
# Eventually overlay stars 
#   Q: Will there be enough there?
# Eventually repoint based on stars
#  ** Before I allow repointing, I should search a star catalog and plot them.


# Now that data are read in, we want to make a 

# Make a plot vs. time showing distance to the ring (i.e., resolution),
# of ring pixels imaged, and 
plt.rcParams['figure.figsize'] = 12, 10 # Make plot normal size

# Convert some things to numpy arrays. Is there any disadvantage to this?

met        = np.array(fits_met)
jd         = np.array(fits_spcutcjd, dtype='d') # 'f' was rounding to one decimal place...
naxis1     = np.array(fits_naxis1)
naxis2     = np.array(fits_naxis2)
target     = np.array(fits_target) # np.array can use string arrays as easily as float arrays
instrument = np.array(fits_spcinst0)
dx_targ    = np.array(fits_spctscx)
dy_targ    = np.array(fits_spctscy)
dz_targ    = np.array(fits_spctscz)
desc       = np.array(fits_reqdesc)
met0       = np.array(fits_startmet)
met1       = np.array(fits_stopmet)
exptime    = np.array(fits_exptime)
rotation   = np.array(fits_spctnaz)
rotation   = np.rint(rotation).astype(int)  # Turn rotation into integer. I only want this to be 0, 90, 180, 270... 
files_short = np.zeros(num_obs, dtype = 'S30')

dist_targ = np.sqrt(dx_targ**2 + dy_targ**2 + dz_targ**2)

phase = np.zeros(num_obs)
utc = np.zeros(num_obs, dtype = 'S30')
et = np.zeros(num_obs)
subsclat = np.zeros(num_obs) # Sub-sc latitude
subsclon = np.zeros(num_obs) # Sub-sc longitude

name_observer = 'New Horizons'
frame = 'J2000'
abcorr = 'LT'

# Fix the MET. The 'MET' field in fits header is actually not the midtime, but the time of the first packet.
# I am going to replace it with the midtime.

met = (met0 + met1) / 2.

# Loop over all images

for i in i_obs:

# Get the ET and UTC, from the JD. These are all times *on s/c*, which is what we want

#  et[i] = cspice.utc2et('JD ' + repr(jd[i]))
  et[i]= 0
#  utc[i] = cspice.et2utc(et[i], 'C', 2)
  utc[i] = 'UT 0'

# Calculate Sun-Jupiter-NH phase angle for each image 

#  (st_jup_sc, ltime) = cspice.spkezr('Jupiter', et[i], frame, abcorr, 'New Horizons') #obs, targ
  st_jup_sc = np.array([0,0,0,0,0,0])
  ltime = 0
#  (st_sun_jup, ltime) = cspice.spkezr('Sun', et[i], frame, abcorr, 'Jupiter')
  st_sun_jup = np.array([0,0,0,0,0,0])
  ltime = 0
#  phase[i] = cspice.vsep(st_sun_jup[0:3], st_jup_sc[0:3])
  phase[i] = 0
  files_short[i] = files[i].split('/')[-1]
# Calc sub-sc lon/lat
  
#  (radius,subsclon[i],subsclat[i]) = cspice.reclat(st_jup_sc[0:3])
  (radius, subsclon[i], subsclat[i]) = (0,0,0) 

# Stuff all of these into a Table

t = Table([i_obs, met, utc, et, jd, files, files_short, naxis1, naxis2, target, instrument, dx_targ, dy_targ, dz_targ, desc, 
           met0, met1, exptime, phase, subsclat, subsclon, naxis1, naxis2, rotation], 
           names = ('#', 'MET', 'UTC', 'ET', 'JD', 'Filename', 'Shortname', 'N1', 'N2', 'Target', 'Inst', 'dx', 'dy', 'dz', 'Desc',
                    'MET Start', 'MET End', 'Exptime', 'Phase', 'Sub-SC Lat', 'Sub-SC Lon', 'dx_pix', 'dy_pix', 'Rotation'))

# Define units for a few of the columns
                    
t['Exptime'].unit = 's'
t['Sub-SC Lat'].unit = 'degrees'

# Create a dxyz_targ column, from dx dy dz. Easy!

t['dxyz'] = np.sqrt(t['dx']**2 + t['dy']**2 + t['dz']**2)

# Get a list of all of the observation descriptions
# ** Note that 'High Phase M/monitoring' is a few of these -- not clearly 'ring' obs

list(set(t['Desc'])) # I can do this for one column, but not two, weirdly... lots of errors.
                     # [got that fixed]

# Get a list of unique observation descriptions. The syntax for this is surprisingly arcane.

keys = (['Target', 'Desc'])
t2 = astropy.table.unique(t, keys = keys)
print t2[keys]

# Do the same, but easier method

keys = (['Target', 'Desc'])
(astropy.table.unique(t, keys = keys))[keys]

# Now do some grouping.  First order them by 'Description'. Make a new table with this.

by_desc = (t.group_by('Desc'))['Desc', 'MET', 'Phase', 'Target']

# Within that table, take the mean of the various quantities

by_desc.groups.aggregate(np.mean)

# Within each Desc field: average Phase and MET, and make a plot of it

plt.plot(by_desc['MET'].groups.aggregate(np.mean), by_desc['Phase'].groups.aggregate(np.mean), ls='None', marker='o')

# Now can I sum the number of records that match?
# Do some filtering: list only the 4x4 data before Mar 1 with latitude < 10 degrees.
#
# &  : bitwise operator, works on NumPy arrays
# && : Does not exist in python (only C)
# and: Not bitwise, and not element-by-element.

# Extract row numbers of the matching records

mask = (t['dx_pix'] == 1024) & (t['Phase'] < 0.4)

# Print this table as latex
astropy.io.ascii.write(t['UTC', 'MET', 'Exptime'][mask],format='latex') # Cool -- it really works!

# Make a plot of phase angle
#quit # NB: IF I uncomment this, then I get a kernel crash later on, when I try to navigate an image. Weird.
plt.plot(met, phase*180/3.14)

plt.plot(t['MET'], t['Phase'])

# On top of this, plot one line per observation
for i in i_obs:
    plt.plot([met[i], met[i]], phase[i]*180/3.14 + np.array([exptime[i], -exptime[i]]), color='red')

plt.xlabel('MET')
plt.ylabel('Phase angle [deg]')
plt.title('Jupiter ring, ' + repr(num_obs) + ' x LORRI')
plt.show()

### MAIN LOOP ####

##########    
# Now start the main loop
##########

i = 0   # Index number to print

#tableheader = "#    File                   Type     Exptime SPCIDNTFY _cr  _rect  _sky  _ext   #"

#print tableheader

IS_DONE = False

r2d = 180 / np.pi

range_selected = range(np.size(files))            # By default, set the range to act upon every file

prompt = "(#)number, (sr)set range of files for action, (q)uit, (l)ist, (n)ext, (p)revious, list (g)roups, \n" + \
         "(h)eader,  (m)edian, (d)s9, (A)rray, (o)utput png, (N)avigate, (" + repr(i) + ") ?: "

t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'].pprint(max_lines=-1, max_width=-1)

columns_print =   ['#', 'MET', 'UTC', 'Exptime', 'Target', 'Rotation', 'Desc']

DO_PLOT_I = False

DO_GROUP = False
    
while (IS_DONE == False): 

    inp = raw_input("File " + repr(i) + ":" + prompt)

    if (inp == 'sr'):
      range_selected = get_range(max = size(files)) # Ask the user to give a range
  
    if (inp == 'q'): # Quit
        IS_DONE = True
        DO_PLOT_I = False
    
    if (inp == 'd'): # Load in DS9, or at least put the command on the screen
        print 'o ' + t['Filename'][i]
        DO_PLOT_I = False

##########
# Display (h)eader
##########
    
    if (inp == 'h'): # Display the header
        f = t['Filename'][i]   
        hdulist = fits.open(f)
        header = hdulist[0].header
        inp2 = raw_input("Search string? ") # If desired, search for a particular string in the header (case-insensitive)
        keys = header.keys()
        if (inp2 == ''):
            print repr(header)
        else:
            for key in keys:
                line = key + ' = ' + repr(header[key])
                if (inp2.upper() in line.upper()):
                    print line
        hdulist.close()
        
##########
# Show (A)rray of images
##########
    
    if (inp == 'A'): # Show array of images. It is not output, but can be cut & pasted from console to textedit.
        range_selected = get_range(maxrange = np.size(files)) # Ask the user to give a range
        num_cols = float(3)
        num_rows = int(math.ceil( np.size(range_selected) / num_cols )) # Round UP
        
        plt.rc('figure', figsize=(15, 15*num_rows/num_cols))
        fig = plt.figure(1)
        plt.rc('figure', figsize=(15, 15*num_rows/num_cols))

        plotnum = 1

        for j in range_selected:
            arr = get_image_nh(file = files[j], polyfit=True)
#            rownum = plotnum / num_cols + 1
#            colnum = plotnum % num_cols + 1        
            plt.subplot(num_rows, num_cols, plotnum)
            plt.imshow(arr)
            plt.title(repr(j) + ', ' + t['Shortname'][j])

            plotnum += 1
            
        plt.show()

##########
# (O)utput multiple images as PNG
##########
            
    if (inp == 'o'): # If a single image, output it. If a range, output them, but all at the same scaling
        range_selected = get_range(maxrange = np.size(files)) # Ask the user to give a range
    
        hdulist = fits.open(t['Filename'][i])
        image = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'
        if hdulist[0].header['TARGET'] == 'IO':
            power = 0               # Only do the polynomial fit for ring (Target = Jupiter)
        else:
            power = 10    
        frac_max = 0.9              # Clip any pixels brighter than this fractional amount  
        mode_0 = scipy.stats.mode(np.round(image.data),axis=None)[0][0] # Get the most popular DN value
        arr_filtered = hbt.remove_brightest(image.data, frac_max)   
        polyfit = hbt.sfit(arr_filtered, power)
        arr_plot = arr_filtered - polyfit
        plt.imshow(arr_plot)
        plt.title(repr(i) + ', ' + t['Shortname'][i])
        plt.show()
        vmin = np.amin(arr_plot)    # Grab the min and max values of these so we can scale others the same way
        vmax = np.amax(arr_plot)
        hdulist.close 
        med_j = np.zeros(np.size(range_selected))
              
        hbt.imsize((5,5))
        for jj,j in enumerate(range_selected):
            hdulist = fits.open(t['Filename'][j])
            image = hdulist['PRIMARY']
            arr = image.data
            mode_j = scipy.stats.mode(np.round(arr), axis=None)[0][0]
            arr *= mode_0 / mode_j # Normalize the entire array, based on the statistical mode

            DO_HISTOGRAM = False
            if DO_HISTOGRAM:            
                bins = np.arange(100,600,1)
                h = np.histogram(arr,bins=bins)
                plt.hist(h,bins=bins)
                plt.plot(bins[0:-1],h[0])
                plt.xlim((220,650))
                plt.title(repr(j) + ', ' + t['Shortname'][j])
                plt.show()
            
            med_j[jj] = np.median(arr)
            plt.imshow(arr - polyfit, vmin=vmin, vmax=vmax)
            plt.title(repr(j) + ', ' + t['Shortname'][j])
            plt.show()
            hdulist.close()

##########
# Enable (m)edian subtraction
##########

    if (inp == 'm'): # Enable median subtraction
        inp2 = raw_input("Enable median subraction?")
        if (inp2 == 'y'):
            DO_MEDIAN_SUBTRACT = 1
        else:
            DO_MEDIAN_SUBTRACT = 0
        
        if (DO_MEDIAN_SUBTRACT):
            range_selected = get_range(maxrange = np.size(files)) # Ask the user to give a range
            num_images = np.size(range_selected)
            image_arr = np.zeros((num_images, 1024, 1024))
            mode_arr = np.zeros(num_images)
            mode = np.zeros(num_images)
            
            for ii,j in enumerate(range_selected):    
                hdulist = fits.open(files[j])
                header = hdulist[0].header 
                image_j = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'
                mode_arr[ii] = scipy.stats.mode(np.round(image_j.data), axis=None)[0][0]
                r = match_image_brightness(image_arr[0], image_j.data, frac=0.5)
                image_arr[ii,:,:] = image_j.data * r[0] + r[1]
                hdulist.close()
                
            med = np.median(image_arr,0) # Median works ok, but if the ring stays in same area, median *is* the ring
            max = np.amax(image_arr,0)   # Max can be contaminated by star light
            min = np.amin(image_arr,0)   # Max can be contaminated by star light
            for ii in range(num_images):
                mode[ii] = scipy.stats.mode(np.round(image_arr[i,:,:]),axis=None)[0][0]
            
#            mode = scipy.stats.mode(np.round(imagearr, axis=1)
            
            plt.imshow(med)
            plt.title('Median, images ' + repr(med_start) + ' .. ' + repr(med_end))
            plt.show()

            plt.imshow(min)
            plt.title('Min, images ' + repr(med_start) + ' .. ' + repr(med_end))
            plt.show()

            plt.imshow(max)
            plt.title('Max, images ' + repr(med_start) + ' .. ' + repr(med_end))
            plt.show()            
            
            plt.imshow(imagearr[0] - med)
            plt.title('Image ' + repr(i) + ' - median')
            plt.show()

##########
# Starting using a (g)roup
##########
        
    if (inp == 'g'): # Start using a Group (ie, )
        t_by_desc = t.group_by('Desc')
        descs = np.array(t_by_desc.groups.keys)
        print '*** Groups ***'
        for i,d in enumerate(descs):
            print repr(i) + ' ' + d[0]
        print "# per group:" + repr(np.diff(t_by_desc.groups.indices))
        inp2 = raw_input("Group # (0 .. " + repr(np.size(t_by_desc.groups.keys)-1) + ")")

        if (inp2 == ''):   # If we hit <cr>, then return to the full list of files
            DO_GROUP = False
            inp = 'l'
            
        if hbt.is_number(inp2):
            num_group = int(inp2)
            name_group = t_by_desc.groups.keys[num_group][0] # Extract the group name
            indices = t['Desc'] == name_group
            tg = t[indices][columns_print] # Extract members of the group itself. I suspect this is not the easiest way to do this??
            print '*** ' + name_group + ' ***'
            tg[columns_print].pprint(max_width=-1, max_lines=-1)
            inp = repr(tg['#'][0])     # Seed the input string so as to display first member of this group
            DO_PLOT_I = True               # And plot it (but can kill this if we don't want this behavior)
            DO_GROUP = True
            
    if (hbt.is_number(inp)):
        i = int(inp)
        DO_PLOT_I = True

##########
# (p)revious image
##########
        
    if (inp == 'p'): # Previous
        if DO_GROUP:
            indices = (np.array(tg['#'])) # Get a list of indices in this group
            index = (indices == i)
            i = int(indices[np.roll(index,-1)])
        
        else:
            i -= 1
            
        DO_PLOT_I = True

##########
# (n)ext image
##########
    
    if (inp == 'n'): # Next image
        if DO_GROUP:
            indices = (np.array(tg['#'])) # Get a list of indices in this group
            index = (indices == i)        # Find index of current one
            i = int(indices[np.roll(index,1)])    # And move backwards in the list

        else:
            i += 1
        DO_PLOT_I = True

##########
# (l)ist all frames
##########
       
    if (inp == 'l'): # List all frames
        t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'].pprint(max_lines=-1, max_width=-1)
        t[columns_print].pprint(max_lines=-1, max_width=-1)

#########
# (N)avigate the image -- plot rings and everything on it
#########

    if (inp == 'N'): # Navigate the image -- plot rings and everything on it
        
# Now look up positions of stars in this field, from a star catalog

        w = WCS(t['Filename'][i])                  # Look up the WCS coordinates for this frame

        center  = w.wcs.crval  # degrees
        name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1'
        print 'got 1'
        stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
        print 'got 2'
        ra_stars  = np.array(stars.array['RAJ2000'])*d2r # Convert to radians
        print 'got 3'
        dec_stars = np.array(stars.array['DEJ2000'])*d2r # Convert to radians
        table_stars = Table(stars.array.data)

# Get an array of points along the ring

        ra_ring, dec_ring = get_jring_points_radec(et)                    # Return as radians              
        x_ring, y_ring    = w.wcs_world2pix(ra_ring*r2d, dec_ring*r2d, 0) # Convert to pixels
        
# Look up velocity of NH, for stellar aberration
        
        et = t['ET'][i]
        abcorr = 'LT+S'
        frame = 'J2000'
#        st,ltime = cspice.spkezr('New Horizons', et, frame, abcorr, 'Sun') # Get velocity of NH 
        st = np.array([0,0,0,0,0,0])
        ltime = 0
        vel_sun_nh_j2k = st[3:6]
        
# Correct stellar RA/Dec for stellar aberration

        radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
        radec_stars_abcorr = hbt.correct_stellab(radec_stars, vel_sun_nh_j2k) # Store as radians

# Convert ring RA/Dec for stellar aberration

        radec_ring        = np.transpose(np.array((ra_ring,dec_ring)))
        radec_ring_abcorr = hbt.correct_stellab(radec_ring, vel_sun_nh_j2k) # radians
        
# Convert RA/Dec values back into pixels
        
        x_stars,        y_stars          = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)  # final arg: row vs column order      
        x_stars_abcorr, y_stars_abcorr   = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
        x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)

        points_stars        = np.transpose((x_stars, y_stars))
        points_stars_abcorr = np.transpose((x_stars_abcorr, y_stars_abcorr))

# Read the image file from disk

        image_polyfit = get_image_nh(t['Filename'][i], frac=1., polyfit = True)
        image_raw     = get_image_nh(t['Filename'][i], raw = True)

# Use DAOphot to search the image for stars. It works really well.

        points_phot = find_stars(image_polyfit)
        
# Now look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
# I call this 'opnav'. It is returned in order (y,x) because that is what imreg_dft uses -- even though it is a bit weird.

        (dy_opnav, dx_opnav) = calc_offset_points(points_phot, points_stars_abcorr, np.shape(image_raw), plot=True)

# Now convert this pixel offset to a radec offset, and tweak wcs.
        
# Now assemble it all into a single composite image
# Remove most of the border -- seee http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
#  ** First plot
#  Plot the image
#  Plot the photometric stars, in boxes
#
#  ** Second plot
#  Plot the catalog stars *plus derived offset* in circles

        xrange = (0, np.shape(image_raw)[0])
        yrange = (0, np.shape(image_raw)[1])

        figs = plt.figure()
        ax1 = figs.add_subplot(1,2,1) # nrows, ncols, plotnum. Returns an 'axis'

        fig1 = plt.imshow(np.log(image_polyfit))

        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        
# Overlay photometric stars
        plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', markersize=12)

# Now make second plot    
        ax2 = figs.add_subplot(1,2,2)
        fig2 = plt.imshow(np.log(image_polyfit))

# Overlay photometric stars
        plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', color='red', markersize=12)

# Overlay catalog stars
#        plt.plot(points_stars[:,0],  points_stars[:,1],  marker='o', ls='None', color='lightgreen')

# Overlay catalog stars, corrected for opnav offset
        plt.plot(points_stars_abcorr[:,0] + dx_opnav, points_stars_abcorr[:,1] + dy_opnav, marker='o', ls='None', color='lightgreen')
        
# Overlay rings points
        plt.plot(x_ring, y_ring, marker = 'o', color='red', ls='-')
 
# Now need to apply stellar aberration (and LT?) to the ring points        
       
        plt.xlim(xrange)
        plt.ylim(yrange)
        
#        fig2.axes.get_xaxis().set_visible(False) # Turn off the axes labels if requested
#        fig2.axes.get_yaxis().set_visible(False)         
#        plt.axis('off') # Suppress all axis, labels, etc. 
#        plt.Axes(figs, [0,0,1,1]) # Axes gets passed a figure, not an axis.
        plt.show()     

# Now assemble it all into a single composite image
# Remove most of the border -- seee http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection=w) # Project using the current WCS projection
#        plt.axis('off') # Suppress all axis, labels, etc. 
        ax = plt.Axes(fig, [0,0,1,1]) 
        fig2 = plt.imshow(np.log(image_polyfit))
        
        plt.plot(points_stars[:,0] + dx_opnav, points_stars[:,1] + dy_opnav, marker='o', ls='None', color='lightgreen', ms=12, mew=1)

        plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', color='pink')

# Plot the ring. For fun, plot it twice, showing effect of stellar aberration

        plt.plot(x_ring + dx_opnav, y_ring + dy_opnav, marker = 'o', color='red', ls='-', ms=8)
        plt.plot(x_ring_abcorr + dx_opnav, y_ring_abcorr + dy_opnav, marker = 'o', color='red', ls='-', ms=8)

        fig2.axes.get_xaxis().set_visible(False)
        fig2.axes.get_yaxis().set_visible(False) 
        plt.xlim((0,1000))
        plt.ylim((0,1000))
        
#        x, y = w.wcs_world2pix(ra, dec, 0)
#        plot(x, y, marker='o', ls='None')        
        plt.show()

######## TEST WCS OFFSETTING ##########

        do_test_wcs = False
        
        if (do_test_wcs):
            w = WCS(t['Filename'][i])                  # Look up the WCS coordinates for this frame
    #        rcParams['figure.figsize'] = 20, 10
            plt.rc('figure', figsize=(16,8))
            figs = plt.figure()
            ax1 = figs.add_subplot(1,2,1) # nrows, ncols, plotnum. Returns an 'axis'
            fig1 = plt.imshow(np.log(image_polyfit))
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)  # final arg: row vs column order      
            plt.plot(x_stars, y_stars, marker='o', ls='None', color='lightgreen', ms=12, mew=1)
                    
            plt.xlim((0,1000))
            plt.ylim((0,1000))
            plt.title('WCS, center = ' + repr(w.wcs.crval))
            plt.xlabel('Dec, center = ' + repr(w.wcs.crval[1]))
            plt.ylabel('RA, center = ' + repr(w.wcs.crval[0]))
       
# Now convert the dx / dy offset (from opnav) into an RA/Dec offset.
# Put the new correct center position into WCS, so plots from now on will be correct.

        do_offset_wcs = False
        
        if (do_offset_wcs):
            m = w.wcs.piximg_matrix
            w.wcs.crval[0] -= (dy_opnav * (m[0,0] + m[1,0])) # RA. I am sort of guessing about these params, but looks good.
            w.wcs.crval[1] -= (dx_opnav * (m[1,1] + m[0,1])) # Dec
    
            ax1 = figs.add_subplot(1,2,2, projection=w) # nrows, ncols, plotnum. Returns an 'axis'
            fig1 = plt.imshow(np.log(image_polyfit))
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)  # final arg: row vs column order      
            
            plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', markersize=12)
    
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)  # final arg: row vs column order      
            x_stars_abcorr, y_stars_abcorr = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
            x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)
    
    
            plt.plot(x_stars, y_stars, marker='o', ls='None', color='lightgreen', mfc = 'None', mec = 'red', ms=12, mew=2, label='Cat')
            plt.plot(x_stars_abcorr, y_stars_abcorr, marker='o', ls='None', mfc = 'None', mec = 'green', ms=12, mew=2, label = 'Cat, abcorr')
            plt.plot(x_ring_abcorr, y_ring_abcorr, marker='o', ls='-', color='blue', mfc = 'None', mec = 'blue', ms=12, mew=2, 
                     label = 'Ring, LT+S')
                    
            plt.xlim((0,1000))
            plt.ylim((0,1000))
            plt.title('WCS, center = ' + repr(w.wcs.crval))
            
        DO_PLOT_I = False
    
##########    
# Load and plot the image, if appropriate flag is set.
# This runs at end of loop.
##########
    
    if (DO_PLOT_I):     
        print t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'][i]

        arr = get_image_nh(t['Filename'][i], polyfit=True)
#        plt.imshow(np.log(arr))
        plt.imshow(arr)
        plt.title(t['Shortname'][i])
        plt.show()
        
#        plt.show()

# Plot Io's position, as a test

#        (state_io, ltime) = cspice.spkezr('Io', et[i], 'J2000', 'LT+S', 'New Horizons')
#        (dist_io, ra_io, dec_io) = cspice.recrad(state_io[0:3])
#        (x_io, y_io) = w.wcs_world2pix(ra_io*r2d, dec_io*r2d, 0)
#        plt.plot(x_io, y_io, mec = 'purple', mew = 4, ms=30, marker = 'o', label='Io', mfc='none')
#        plt.legend()
#        plt.show()
        

#    
# Define a method to list all of the images, or a sub-range of them.

# class NH_Jring:
# Methods: get images
# populate fields
# add new fields
# compute geometry
# compute backplane for an image
# display a single image (in window, or inline)
# Overlay a ring
#
# Q: Is the object *one* image? Or a set of 571?
# RingImages->Load(list)
# RingImages->utc
# RingImages->met
#
# Or should we just use a structure
#
# Other functionality;
#  Select an image or set of images, either by clicking individually and/or filtering by criteria
#  Process the selected image(s)
#  View the selected image(s)

quit  # When the IS_DONE loop is finished, execution drops to here

subsollon  = np.empty_like(jd)

name_observer = 'Earth'
name_target   = 'Pluto'
method = 'Near Point'
fixref = 'IAU_PLUTO'
abcorr = 'LT'
r2d    = 360. / 2. / math.pi


def calc_offset_slow(image_phot):
    "Find the offset between the images."
    "This uses a painfully slow brute force approach. Definitely not recommended."
    "Hasn't been tested as a function"
    
    corr = np.zeros((60,60))
    for i in range(60):
        for j in range(60):
            corr[i,j] = np.sum(np.maximum(im1, np.roll(np.roll(im0,i-30,axis=0),j-30,axis=1)))

# Normalize the correlation array so there is one peak, with range 0 .. 1
     
    corr = corr - np.min(corr)
    corr = np.max(corr) - corr
    corr = corr / np.max(corr)

# Use DAOfind to locate this peak. Could do it other methods, but this works well, even if overkill.
    
    s2 = daofind(corr, fwhm=10, threshold=0.8)

# Now get the offsets. This is in pixels.

    dx = int(np.round(s2['xcentroid'][0]))
    dy = int(np.round(s2['ycentroid'][0]))
    
# And roll the star catalog image

    im0_rolled = np.roll(np.roll(image_0,dy-30,axis=0), dx-30, axis=1)        
    
    return im0_rolled
    
