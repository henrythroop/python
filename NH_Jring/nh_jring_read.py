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
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   photutils import daofind
import wcsaxes
import HBT as hbt
import warnings
import imreg_dft as ird

#function calc_median_framelist(t)


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

#        ax3 = figs.add_subplot(2,2,3)
        
        plt.show()
        
# Plot the pair of images

#        plt.imshow(image_phot + 
#          np.roll(np.roll(image_from_list_points(points_cat, np.shape(image.data),diam_kernel),t0[0],0),t0[1],1))
#        plt.show()
        
    return t0
        
hbt.set_plot_defaults()
  
d2r = np.pi /180.
r2d = 1. / d2r
         
dir_data = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'
# Start up SPICE

file_tm  = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"
cspice.furnsh(file_tm)

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
rotation   = np.rint(rotation).astype(int)  # Turn rotation into integer. I only want this to be 0, 90, 180, 270... I don't care about the resolution so much.

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

  et[i] = cspice.utc2et('JD ' + repr(jd[i]))
  utc[i] = cspice.et2utc(et[i], 'C', 2)

# Calculate Sun-Jupiter-NH phase angle for each image 

  (st_jup_sc, ltime) = cspice.spkezr('Jupiter', et[i], frame, abcorr, 'New Horizons') #obs, targ
  (st_sun_jup, ltime) = cspice.spkezr('Sun', et[i], frame, abcorr, 'Jupiter')
  phase[i] = cspice.vsep(st_sun_jup[0:3], st_jup_sc[0:3])
  
# Calc sub-sc lon/lat
  
  (radius,subsclon[i],subsclat[i]) = cspice.reclat(st_jup_sc[0:3])

# Stuff all of these into a Table

t = Table([i_obs, met, utc, et, jd, files, naxis1, naxis2, target, instrument, dx_targ, dy_targ, dz_targ, desc, 
           met0, met1, exptime, phase, subsclat, subsclon, naxis1, naxis2, rotation], 
           names = ('#', 'MET', 'UTC', 'ET', 'JD', 'Filename', 'N1', 'N2', 'Target', 'Inst', 'dx', 'dy', 'dz', 'Desc',
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
quit
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
         "(h)eader,  (m)edian, (d)s9, (" + repr(i) + ") ?: "

t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'].pprint(max_lines=-1, max_width=-1)

columns_print =   ['#', 'MET', 'UTC', 'Exptime', 'Target', 'Rotation', 'Desc']

DO_PLOT_I = False

DO_GROUP = False
    
while (IS_DONE == False): 

    inp = raw_input("File " + repr(i) + ":" + prompt)

    if (inp == 'sr'):
        inp2 = raw_input("Range of files [e.g. *; 1-10; 44]: ")

        if (inp2.find('-') >= 0):        # Range of files
            min = int((inp2.split('-'))[0])
            max = int((inp2.split('-'))[1])
            range_selected = range(min, max+1)

        elif (inp2.find('*') >= 0):        # All files
            range_selected = range(size(files))

        else:        # Only a single file
            range_selected = [int(inp2)]
            print 'Range: ' + repr(range_selected)
  
    if (inp == 'q'): # Quit
        IS_DONE = True
    
    if (inp == 'd'): # Load in DS9, or at least put the command on the screen
        print 'o ' + t['Filename'][i]
    
    if (inp == 'h'): # Display the header
#        f = t['Filename'][i]   
#        hdulist = fits.open(f)
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
    
    if (inp == 'f'):
        pass
    
    if (inp == 'm'): # Enable median subtraction
        inp2 = raw_input("Enable median subraction?")
        if (inp2 == 'y'):
            DO_MEDIAN_SUBTRACT = 1
        else:
            DO_MEDIAN_SUBTRACT = 0
        
        if (DO_MEDIAN_SUBTRACT):
            med_start = int(raw_input("First frame as median: "))
            med_end   = int(raw_input("Last frame as median: "))
            imagearr = np.zeros((med_end - med_start, 1024, 1024))
            for j in np.arange(med_start, med_end):    
                hdulist = fits.open(files[j])
                header = hdulist[0].header 
                imagej = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'
                imagearr[j-med_start,:,:] = imagej.data
            med = np.median(imagearr,0)
            plt.imshow(med)
            plt.title('Median, images ' + repr(med_start) + ' .. ' + repr(med_end))
            plt.show()
            
            plt.imshow(image.data - med)
            plt.title('Image ' + repr(i) + ' - median')
            plt.show()
        
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
        
    if (inp == 'p'): # Previous
        if DO_GROUP:
            indices = (np.array(tg['#'])) # Get a list of indices in this group
            index = (indices == i)
            i = int(indices[np.roll(index,-1)])
        
        else:
            i += 1
            
        DO_PLOT_I = True
    
    if (inp == 'n'):
        if DO_GROUP:
            indices = (np.array(tg['#'])) # Get a list of indices in this group
            index = (indices == i)        # Find index of current one
            i = int(indices[np.roll(index,1)])    # And move backwards in the list

        else:
            i -= 1
        DO_PLOT_I = True
        
    if (inp == 'l'):
        t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'].pprint(max_lines=-1, max_width=-1)
        t[columns_print].pprint(max_lines=-1, max_width=-1)
    
# If the flag is set, then load and plot the image
    
    if (DO_PLOT_I):     
#        t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'][i].pprint(max_lines=-1, max_width=-1)
        print t['#', 'MET', 'UTC', 'Exptime', 'Target', 'Desc'][i]

        if ('hdulist') in locals():             # If there is already an hdulist, then close it. (Might double-close it, but that is OK.)
            hdulist.close()
        f = t['Filename'][i]   
        hdulist = fits.open(f)
        image = hdulist['PRIMARY'] # options are 'PRIMARY', 'LORRI Error', 'LORRI Quality'

# Fit the data using a polynomial

        if t['Target'][i] == 'IO':
            power = 0               # Only do the polynomial fit for ring (Target = Jupiter)
        else:
            power = 10
       
        p = hbt.sfit(image.data, power)

# Now calculate the ring points...

        radii_ring = np.array([1220000., 129000.])  # Inner and outer radius to plot
        num_pts_ring = 100
        
        ring_lon = np.linspace(0, 2. * np.pi, num_pts_ring)
        ra_ring  = np.zeros(num_pts_ring)
        dec_ring = np.zeros(num_pts_ring)
        
        frame = 'J2000'
        abcorr = 'LT'
        rot = cspice.pxform('IAU_Jupiter', frame, t['ET'][i]) # Get matrix from arg1 to arg2
        
        st,ltime = cspice.spkezr('Jupiter', t['ET'][i], frame, abcorr, 'New Horizons')
        pos = st[0:3]
        vel = st[3:6] # velocity, km/sec, of jupiter
        
        st,ltime = cspice.spkezr('New Horizons', t['ET'][i], frame, abcorr, 'Sun') # Get velocity of NH 
        vel_sun_nh_j2k = st[3:6]
        
        quit      
        for radius_ring in radii_ring:
            for j in range(num_pts_ring):
                xyz = np.zeros(3)
                xyz = np.array((radius_ring * np.cos(ring_lon[j]), radius_ring * np.sin(ring_lon[j]), 0.))
                
                d_j2000_xyz = np.dot(rot,xyz)  # Manual says that this is indeed matrix multiplication
                j2000_xyz = 0 * d_j2000_xyz
                j2000_xyz[0:3] = d_j2000_xyz[0:3] # Looks like this is just copying it

                rho_planet    = pos                     # Position of planet
                rho_ring      = rho_planet + j2000_xyz  # Vector obs-ring
                dist_ring     = cspice.vnorm(rho_ring)*1000 # Convert to km... CHECK UNITS!
                
                range_out, ra, dec = cspice.recrad(rho_ring) # 'range' is a protected keyword in python!
                
                ra_ring[j] = ra     # save RA, Dec as radians
                dec_ring[j] = dec
                
        w = WCS(t['Filename'][i])                  # Look up the WCS coordinates for this frame

        x_ring, y_ring = w.wcs_world2pix(ra_ring*r2d, dec_ring*r2d, 0)
        
# Now look up positions from a star catalog

        center  = w.wcs.crval  # degrees
        name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1'
        stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=True, catalog_db = name_cat)
        ra_cat = np.array(stars.array['_RAJ2000'])*d2r # Convert to radians
        dec_cat = np.array(stars.array['_DEJ2000'])*d2r # Convert to radians
        table_stars = Table(stars.array.data)
        
        DO_PLOT_I = False

# Correct stellar RA/Dec for stellar aberration

        radec_cat = np.transpose(np.array((ra_cat,dec_cat)))
        radec_cat_abcorr = hbt.correct_stellab(radec_cat, vel_sun_nh_j2k) # Store as radians

# Convert ring RA/Dec for stellar aberration

        radec_ring = np.transpose(np.array((ra_ring,dec_ring)))
        radec_ring_abcorr = hbt.correct_stellab(radec_ring, vel_sun_nh_j2k) # radians
        
# Convert RA/Dec values back into pixels
        
        x_cat,        y_cat        = w.wcs_world2pix(radec_cat[:,0]*r2d,   radec_cat[:,1]*r2d, 0)  # final arg deal with row vs column order      
        x_cat_abcorr, y_cat_abcorr = w.wcs_world2pix(radec_cat_abcorr[:,0]*r2d, radec_cat_abcorr[:,1]*r2d, 0)
        x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)

        points_cat = np.transpose((x_cat, y_cat))
        points_cat_abcorr = np.transpose((x_cat_abcorr, y_cat_abcorr))

# Use DAOphot to search the image for stars. It works really well.
        
        points_phot = find_stars(image.data - p)

# Now look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
# I call this 'opnav'. It is returned in order (y,x) because that is what imreg_dft uses -- even though it is a bit weird.

        (dy_opnav, dx_opnav) = calc_offset_points(points_phot, points_cat_abcorr, np.shape(image.data), plot=True)

# Now convert this pixel offset to a radec offset, and tweak wcs.

        quit
        
# Now assemble it all into a single composite image
# Remove most of the border -- seee http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces
#  ** First plot
#  Plot the image
#  Plot the photometric stars, in boxes
#
#  ** Second plot
#  Plot the catalog stars *plus derived offset* in circles

        xrange = (0, np.shape(image.data)[0])
        yrange = (0, np.shape(image.data)[1])

        figs = plt.figure()
        ax1 = figs.add_subplot(1,2,1) # nrows, ncols, plotnum. Returns an 'axis'

        fig1 = plt.imshow(np.log(image.data - p))

        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        
# Overlay photometric stars
        plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', markersize=12)

# Now make second plot    
        ax2 = figs.add_subplot(1,2,2)
        fig2 = plt.imshow(np.log(image.data - p))

# Overlay photometric stars
        plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', color='red', markersize=12)

# Overlay catalog stars
#        plt.plot(points_cat[:,0],  points_cat[:,1],  marker='o', ls='None', color='lightgreen')

# Overlay catalog stars, corrected for opnav offset
        plt.plot(points_cat_abcorr[:,0] + dx_opnav, points_cat_abcorr[:,1] + dy_opnav, marker='o', ls='None', color='lightgreen')
        
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
        fig2 = plt.imshow(np.log(image.data - p))
        
        plt.plot(points_cat[:,0] + dx_opnav, points_cat[:,1] + dy_opnav, marker='o', ls='None', color='lightgreen', ms=12, mew=1)

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
            fig1 = plt.imshow(np.log(image.data - p))
            x_cat,        y_cat        = w.wcs_world2pix(radec_cat[:,0]*r2d,   radec_cat[:,1]*r2d, 0)  # final arg deal with row vs column order      
            plt.plot(x_cat, y_cat, marker='o', ls='None', color='lightgreen', ms=12, mew=1)
                    
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
            fig1 = plt.imshow(np.log(image.data - p))
            x_cat,        y_cat        = w.wcs_world2pix(radec_cat[:,0]*r2d,   radec_cat[:,1]*r2d, 0)  # final arg deal with row vs column order      
            
            plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', markersize=12)
    
            x_cat,        y_cat        = w.wcs_world2pix(radec_cat[:,0]*r2d,   radec_cat[:,1]*r2d, 0)  # final arg deal with row vs column order      
            x_cat_abcorr, y_cat_abcorr = w.wcs_world2pix(radec_cat_abcorr[:,0]*r2d, radec_cat_abcorr[:,1]*r2d, 0)
            x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)
    
    
            plt.plot(x_cat, y_cat, marker='o', ls='None', color='lightgreen', mfc = 'None', mec = 'red', ms=12, mew=2, label='Cat')
            plt.plot(x_cat_abcorr, y_cat_abcorr, marker='o', ls='None', mfc = 'None', mec = 'green', ms=12, mew=2, label = 'Cat, abcorr')
            plt.plot(x_ring_abcorr, y_ring_abcorr, marker='o', ls='-', color='blue', mfc = 'None', mec = 'blue', ms=12, mew=2, 
                     label = 'Ring, LT+S')
                    
            plt.xlim((0,1000))
            plt.ylim((0,1000))
            plt.title('WCS, center = ' + repr(w.wcs.crval))
        
#        plt.show()

# Plot Io's position, as a test

#        (state_io, ltime) = cspice.spkezr('Io', et[i], 'J2000', 'LT+S', 'New Horizons')
#        (dist_io, ra_io, dec_io) = cspice.recrad(state_io[0:3])
#        (x_io, y_io) = w.wcs_world2pix(ra_io*r2d, dec_io*r2d, 0)
#        plt.plot(x_io, y_io, mec = 'purple', mew = 4, ms=30, marker = 'o', label='Io', mfc='none')
#        plt.legend()
#        plt.show()
        
# Draw the images

# Now create an image (off-screen) of the catalog stars only, and the DAOphot stars only

#        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
#        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))    
        quit


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

quit

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
    
