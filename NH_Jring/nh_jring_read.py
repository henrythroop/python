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
import os.path
from   subprocess import call
import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
import numpy as np
import astropy.modeling
from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   skimage.transform import resize
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   photutils import daofind
import HBT as hbt

#function calc_median_framelist(t)

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
num_obs = size(files)
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
rcParams['figure.figsize'] = 12, 10 # Make plot normal size

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

dist_targ = sqrt(dx_targ**2 + dy_targ**2 + dz_targ**2)

phase = np.zeros(num_obs)
utc = np.zeros(num_obs, dtype = 'S30')
et = np.zeros(num_obs)
subsclat = np.zeros(num_obs) # Sub-sc latitude
subsclon = np.zeros(num_obs) # Sub-sc longitude

name_observer = 'New Horizons'
frame = 'J2000'
abcorr = 'LT'

# Stuff all of this data 
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
#  sbusclon[i] = lon
#  subsclat[i] = lat

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

range_selected = range(size(files))            # By default, set the range to act upon every file

prompt = "(#)number, (sr)set range of files for action, (q)uit, (l)ist, (n)ext, (p)revious, list (g)roups, \n" + \
         "(h)eader,  (d)s9, (" + repr(i) + ") ?: "

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
        
    if (inp == 'g'): # Start using a Group (ie, )
        t_by_desc = t.group_by('Desc')
        descs = np.array(t_by_desc.groups.keys)
        print '*** Groups ***'
        for i,d in enumerate(descs):
            print repr(i) + ' ' + d[0]
        print "# per group:" + repr(np.diff(t_by_desc.groups.indices))
        inp2 = raw_input("Group # (0 .. " + repr(size(t_by_desc.groups.keys)-1) + ")")

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
            i = int(indices[roll(index,-1)])
        
        else:
            i += 1
            
        DO_PLOT_I = True
    
    if (inp == 'n'):
        if DO_GROUP:
            indices = (np.array(tg['#'])) # Get a list of indices in this group
            index = (indices == i)        # Find index of current one
            i = int(indices[roll(index,1)])    # And move backwards in the list

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
#            p = fit_p(p_init, x, y, image.data)
        
#        imshow(skimage.exposure.equalize_hist(image.data - p), cmap=get_cmap('Greys'))

#        plt.show()

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
                
                ra_ring[j] = ra
                dec_ring[j] = dec
                
        w = WCS(t['Filename'][i])                  # Look up the WCS coordinates for this frame

        x_ring, y_ring = w.wcs_world2pix(ra_ring*r2d, dec_ring*r2d, 0)
        
# Now look up positions from a star catalog

        center  = w.wcs.crval  # degrees
        name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1'
        stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=True, catalog_db = name_cat)
        ra = np.array(stars.array['_RAJ2000'])
        dec = np.array(stars.array['_DEJ2000'])
        table_stars = Table(stars.array.data)
        
        ra_abcorr = ra.copy() # Degrees
        dec_abcorr  = dec.copy()
        
        x_cat, y_cat = w.wcs_world2pix(ra, dec, 0)
        
        DO_PLOT_I = False

# Correct for stellar aberration

        x_cat_abcorr = x_cat.copy()
        y_cat_abcorr = y_cat.copy()

        for i in range(len(x_cat)):
            pos_i = cspice.radrec(1., ra[i]*d2r, dec[i]*d2r)
            pos_i_abcorr = cspice.stelab(pos_i, vel_sun_nh_j2k)
            rang, ra_abcorr[i], dec_abcorr[i] = cspice.recrad(pos_i_abcorr)

        ra_abcorr *= r2d # degrees
        dec_abcorr *= r2d
        x_cat_abcorr, y_cat_abcorr = w.wcs_world2pix(ra_abcorr, dec_abcorr, 0)
        
# Use DAOphot to search the image for stars. It works really well.
        
        mean, median, std = sigma_clipped_stats(image.data - p, sigma=3.0, iters=5)
        sources = daofind(image.data - p - median, fwhm=3.0, threshold=5.*std)
        x_phot = sources['xcentroid']
        y_phot = sources['ycentroid']

# Now assemble it all into a single composite image
# Remove most of the border -- seee http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.set_cmap('Greys')
        plt.axis('off') # Suppress all axis, labels, etc. 
        ax = plt.Axes(fig, [0,0,1,1]) 
        fig2 = plt.imshow(np.log(image.data - p), cmap=get_cmap('Greys'), interpolation='nearest')
        plt.plot(x_phot, y_phot, marker='o', ls='None')
        plt.plot(x_cat, y_cat, marker='o', ls='None', color='lightgreen')
        plt.plot(x_cat_abcorr, y_cat_abcorr, marker='o', ls='None', color='lightgreen')

        plt.plot(x_ring,y_ring, marker = 'o', color='red', ls='-')
        fig2.axes.get_xaxis().set_visible(False)
        fig2.axes.get_yaxis().set_visible(False) 
        plt.xlim((0,1000))
        plt.ylim((0,1000))

        quit
      
# Create new arrays with photometric and VO star images in them.
        xx, yy = np.mgrid[:10, :10]
        circle = (xx - 4.5) ** 2 + (yy - 4.5) ** 2
        kernel = roll(roll(circle,5,axis=0),5,axis=1)
        
        image_phot = np.zeros(shape(image.data))
        image_cat  = np.zeros(shape(image.data))
        for i in range(len(x_cat)):                 # XXX Why do I have to use index variable, not for x, y??
            x = int(x_cat[i])
            y = int(y_cat[i])
            if (x > 0) & (x < 950) & (y > 0) & (y < 950):
                image_cat[x:x+10, y:y+10] = kernel

        for i in range(len(x_phot)):
            x = int(x_phot[i])
            y = int(y_phot[i])
            if (x > 0) & (x < 950) & (y > 0) & (y < 950):
                image_phot[x:x+10, y:y+10] = kernel     
    
        corr = np.zeros((100,100))
        for i in range(100):
            for j in range(100):
                corr[i,j] = np.sum(np.maximum(image_phot, roll(roll(image_cat,i-50,axis=0),j-50,axis=1)))

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

        image_cat_rolled = roll(roll(image_cat,dy-50,axis=0), dx-50, axis=1)        

# Now assemble it all into a single composite image
# Remove most of the border -- seee http://stackoverflow.com/questions/9295026/matplotlib-plots-removing-axis-legends-and-white-spaces

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.set_cmap('Greys')
        plt.axis('off') # Suppress all axis, labels, etc. 
        ax = plt.Axes(fig, [0,0,1,1]) 
        fig2 = plt.imshow(np.log(image.data - p), cmap=get_cmap('Greys'), interpolation='nearest')
        plt.plot(x_cat + (dy-50), y_cat + (dx-50), marker='o', ls='None', color='lightgreen', ms=12, mew=1)
        plt.plot(x_phot, y_phot, marker='o', ls='None', color='pink')

        plt.plot(x_ring + (dy-50)-24,y_ring + (dx-50)-32, marker = 'o', color='red', ls='-', ms=8)
        fig2.axes.get_xaxis().set_visible(False)
        fig2.axes.get_yaxis().set_visible(False) 
        plt.xlim((0,1000))
        plt.ylim((0,1000))
        
#        x, y = w.wcs_world2pix(ra, dec, 0)
#        plot(x, y, marker='o', ls='None')        
        show()

# Now create an image (off-screen) of the catalog stars only, and the DAOphot stars only
        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')

        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))    
        quit

#quit

        fig = plt.figure()
        ax = fig.add_subplot(111)        
        rnd = np.random.random((100,100))
        plt.imshow(rnd, interpolation='nearest')
        ax = plt.Axes(fig, [0,0,1,1]) 
        plt.axis('off') # Suppress all axis, labels, etc. 
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.xlim((0,100))
        plt.ylim((0,100))
        plt.show()
        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        print repr(shape(data))
        fig.savefig('out.png', bbox_inches='tight', pad_inches=0)
        
#        imshow(np.log(image.data - p), cmap=get_cmap('Greys'))
#        plot(x, y, marker='o', ls='None')
#        plt.xlim((0,1000))
#        plt.ylim((0,1000))
#    
            
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



