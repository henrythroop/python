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
import matplotlib.pyplot as plt # pyplot
import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *  # So I can change plot size: 
import cspice
                     # http://stackoverflow.com/questions/332289/how-do-you-change-the-size-of-figures-drawn-with-matplotlib

year_obs = 2015		# Set this to 2014 or 2015. That will indicate which dataset to analyze

if (year_obs == 2015):
  dir_data = "/Users/throop/Dropbox/data/SALT_Pluto_2015/product"

if (year_obs == 2014):
  dir_data = "/Users/throop/Dropbox/data/SALT_Pluto_2014/product"

# Start up SPICE

file_tm  = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"
cspice.furnsh(file_tm)

# Get the full list of files

file_list = glob.glob(dir_data + '/*fits')
files = np.array(file_list)

# Read the JD from each file. Then sort the files based on JD.

jd = []
for file in files:
    hdulist = fits.open(file)
    jd.append(hdulist[0].header['JD'])
    hdulist.close()
    
indices = np.argsort(jd)
files = files[indices] # Resort by JD

fits_exptime = []	# new list (same as array) 
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

fits_lampid  = []

# print hdulist[0].header to display the whole thing
# if OBSMODE = IMAGING or INSTRUME=SALTICAM, then no gratings.

files_short = np.array(files)
for i in range(files.size):
    files_short = files[i].split('/')[-1]  # Get just the filename itself

for file in files:
    print "Reading file " + file

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
       
    hdulist.close() # Close the FITS file

print object
print "done"

rcParams['figure.figsize'] = 12, 10 # Make plot normal size

files_short = []
is_raw_fits = []		# Flag: True if file is a 'raw' fits (no _cr_ext_sky extensions). 
				# False if it's been created in the pipeline.

for file in files:
    files_short.append(file.replace(dir_data + '/', '').replace('.fits',''))
    is_raw_fits.append(len(files_short[-1]) == len('mbxgpS201505030233'))		# True if no extensions
    
# Put this all into a Pandas DataFrame

frame = DataFrame({'file': files_short,
                   'object': fits_object,
                   'exptime': fits_exptime,
                   'date_obs': fits_date_obs,
                   'utc_obs' : fits_utc_obs,
                   'jd' : fits_jd,
                   'obsmode' : fits_obsmode,
                   'moonang' : fits_moonang,
                   'instrument' : fits_instrume,
                   'gainset' : fits_gainset,
                   'rospeed' : fits_rospeed,
                   'ccdsum' : fits_ccdsum,
                   'pixscale' : fits_pixscale,
                   'proposer' : fits_proposer,
                   'propid' : fits_propid,
                   'airmass' : fits_airmass,
                   'ra' : fits_ra,
                   'dec': fits_dec,
                   'telra' : fits_telra,
                   'teldec' : fits_teldec,
                   'grating' : fits_grating,
                   'gr_sta' : fits_gr_sta,
                   'gr_angle' : fits_gr_angle,
                   'masktype' : fits_masktyp,
                   'maskid' : fits_maskid})

# Convert some things to numpy arrays. Is there any disadvantage to this?

jd          = np.array(fits_jd)
airmass     = np.array(fits_airmass)
object      = np.array(fits_object, dtype='S30') # np.array can use string arrays as easily as float arrays
instrument  = np.array(fits_instrume)
proposer    = np.array(fits_proposer)
exptime     = np.array(fits_exptime)
date_obs    = np.array(fits_date_obs)
utc_obs     = np.array(fits_utc_obs)
obsmode     = np.array(fits_obsmode)
lampid      = np.array(fits_lampid)
grating     = np.array(fits_grating)
tilt        = np.array(fits_gr_angle)
tilt_anno   = np.array(fits_gr_angle, dtype = 'S30')		# Annotated tilt, as a string
obsmode     = np.array(fits_obsmode)  # Remove extra spaces

subsollon  = np.empty_like(jd)

name_observer = 'Earth'
name_target   = 'Pluto'
method = 'Near Point'
fixref = 'IAU_PLUTO'
abcorr = 'LT'
r2d    = 360. / 2. / math.pi

# 2014 uses 'HD 146233', while 2015 is 'HD146233'. Probably my error in making Step-2? Make them match.

for i in range(size(object)):
  oi = object[i]
  if oi == "HD146233": object[i] = "HD 146233"

# Calculate Pluto Subsolar longitude for each observation

for i in range(size(jd)):
    et = cspice.utc2et("JD" + repr(jd[i]))
    spoint = cspice.subsol(method, name_target, et, abcorr, name_observer)
    (radius, lon, lat) = cspice.reclat(spoint)
    if (lon < 0):
      lon += 2. * math.pi
    subsollon[i] = lon			# Longitude, in radians

# If Lampid = 'NONE', convert it to '' to make it easier to handle
# NONE is used during target observations
# QTH2 is used during flats

lampid[(lampid == 'NONE') | (lampid == 'QTH2')] = ''

# print it

DataFrame(frame, columns = ['jd', 'file', 'proposer', 'object'])

# Make a plot from the DataFrame

# fig2 = plt.figure()
# ay = fig2.add_subplot()
# plt.plot(frame['jd'].values, frame['ra'].values)

# List all of the observations with Target = Pluto
# Make a plot their JD vs. airmass

# Count up how many images we have of the different types
# To do a uniq(), convert to set and then back to list
# Regular python has an .index() method, but that just returns the first instance. 
# Regular python also has a .count() method.
# For a where() function, use NumPy.

object_uniq = list(set(fits_object))
formatstr = "{0:10}:{1:6}"
print formatstr.format("Object", "#")
print formatstr.format("-----", "---")

for obj in object_uniq:
    print formatstr.format(obj, repr(fits_object.count(obj)))
print

# Make a list of how many different instruments or modes we use

instrument_uniq = list(set(fits_instrume))
print formatstr.format("Instrument", "#")
print formatstr.format("-----", "---")
for inst in instrument_uniq:
    print "{0:10}: {1:6}".format(inst, repr(fits_instrume.count(inst)))
print

# Make a list of all the people whose data I have
proposer_uniq = list(set(fits_proposer))
print formatstr.format("Proposer", "#")
print formatstr.format("-----", "---")
for prop in proposer_uniq:
    print "{0:10}: {1:6}".format(prop, repr(fits_proposer.count(prop)))    
print

# Make a list of all the gratings we use
grating_uniq = list(set(fits_grating))
print formatstr.format("Grating", "#")
print formatstr.format("-----", "---")
for grat in grating_uniq:
    print "{0:10}: {1:6}".format(grat, repr(fits_grating.count(grat)))    
print

# Make a list of all the tilt angles we use
tilt_uniq = sorted(list(set(fits_gr_angle)))
print formatstr.format("Tilt", "#")
print formatstr.format("-----", "---")
for til in tilt_uniq:
    print "{0:10}: {1:6}".format(til, repr(fits_gr_angle.count(til)))    
print

# Make a list of all the dates we have data for, and # of obs on each
date_obs_uniq = sorted(list(set(fits_date_obs)))
print formatstr.format("Date", "#")
print formatstr.format("-----", "---")
for d in date_obs_uniq:
    print "{0:10}: {1:6}".format(d, repr(fits_date_obs.count(d)))    
print    
    
# Make a plot of airmass vs. JD

#rcParams['figure.figsize'] = 12, 30 # Make plot big
#
#plt.plot(jd - min(jd), fits_airmass, 'ro')
#plt.xlabel('JD-' + repr(min(jd)))
#plt.ylabel('Airmass')
#plt.show()

##########

ms = 10 # markersize
mew = 2 # markeredgewidth

# red circle: pluto

rcParams['figure.figsize'] = 12, 30 # Make plot big

plt.plot(jd - min(jd), jd - np.trunc(jd), color='black', marker='o', ms=3, ls='none', label = 'Any')


w = (obsmode == 'FABRY-PEROT')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='white', marker='o', ms=ms+10, mew=mew-1,
         label = 'FABRY-PEROT', ls='none')
         
w = (object == 'Pluto') & (instrument == 'SALTICAM')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='orange', marker='s', 
         mew=mew, ms=ms+2, linestyle='none', label = 'Pluto SALTICAM')

w = (object == 'Pluto') & (instrument == 'RSS')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='red', marker='o', 
         mew=mew, ms=ms+2, linestyle='none', label = 'Pluto RSS')

w  = (object == 'FLAT')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='white', marker='o', 
         markersize=10, linestyle='none', label = 'Flat', ls='none') # white circles: flats

w = ((object == 'HD 146233') & (instrument == 'RSS'))
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='yellow', marker='o', markersize=10, 
         linestyle='none', label = 'HD RSS')# yellow: star

w = ((object == 'HD 146233') & (instrument == 'SALTICAM'))
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='yellow', marker='s', markersize=10, 
         linestyle='none', label = 'HD SALTICAM')
         
w = (object == 'ARC')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), color='green', marker='+', markersize=ms, markeredgewidth=mew, 
         label = 'ARC', ls='none')#   ARC

w = (object == 'BIAS')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), 'b+', markersize=ms, markeredgewidth=mew,
         label = 'BIAS')# BIAS: blue +

         
#w = np.where(instrument == 'SALTICAM')
#plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), 'ws', mew=mew-1, ms=ms-1,
#         label = 'SALTICAM')# SALTICAM: white squares

w = np.where(proposer != 'Throop')
plt.plot(jd[w] - min(jd), jd[w] - np.trunc(jd[w]), 'rx', ms=ms, mew=mew-1,
         label = 'Other PI')# Red X: other ppl's data. Plot this last!
 
plt.xlabel('JD-' + repr(min(np.trunc(jd))),fontsize=15)
plt.ylabel('fp(JD)', fontsize=21)
plt.title('SALT observations, ' + repr(year_obs), fontsize=24)        
         
plt.legend(loc='upper center') # Auto create legend

plt.savefig('salt.png')  # Must save before doing show() -- otherwise output is blank
plt.show()               # force display into iPython console


##########
# Make a plot showing rotational phase vs. time
##########

# Extract just the Pluto spectroscopy observations

rcParams['figure.figsize'] = 8, 6 # Make plot normal size

w = (object == 'Pluto') & (instrument == 'RSS')

plot(jd[w], subsollon[w]*r2d, marker = 'x', color = 'red', linestyle='none')
plt.xlabel('JD')
plt.ylabel('Pluto Sub-Solar Longitude [deg]')
show()

##########
# Read in and plot the actual SALTICAM Pluto images
# These compare properly to ds9 sky survey images: "Invert X; 90 degrees"
# These do not match my finder charts, since those have wrong RA/Dec for Pluto.

rcParams['figure.figsize'] = 12, 10 # Make plot normal size

DO_PLOT_IMAGES = False

w = (np.where(  ((object =='Pluto') | (object =='HD146233')) & (instrument == 'SALTICAM')))[0]
files_w = files[w]

if (DO_PLOT_IMAGES):
    for i in w:
        file = files[i]
        hdulist = fits.open(file)
        image = hdulist['SCI'].data
        plt.imshow(log(image))
        plt.title(file.split('/')[-1] + ', ' + date_obs[i] + ' ' + utc_obs[i][:5] + 
                  ', ' + repr(int(round(exptime[i]))) + ' sec, ' + object[i]) 
        # Just the filename
        plt.set_cmap('hot')  # can also do this p.set_cmap if I would have done p = plt.imshow()
        show()
    
w = (np.where((object == 'Pluto') & (instrument == 'SALTICAM')))[0]

##########
# Annotate the tilt angle with the Red / Blue 
##########

for i in range(size(files)):
  color = ""

  if (tilt[i] == ''):		# Tilt is a string (not a float). Handle case where tilt is invalid. 
    ftilt = 0
  else:
    ftilt = float(tilt[i])

  if (ftilt > 13) & (ftilt < 14):
    color = "BL"
  if (ftilt > 20) & (ftilt < 21):
    color = "R"

  tilt_anno[i] = "{:7} {:2}".format(tilt[i], color)

##########
# Make a list of all of the observations, per day
# To do this, just print by list of int(jd)
##########

file_out = dir_data + "/SALT_Pluto_Observing_Log.txt"

jd_uniq = sorted(np.array(list(set(np.trunc(jd)))))

calfile = "/Users/throop/python/salt_pluto_cal.txt"

night_num = 1

lines = []

# Convert tilt from string to number. '' -> 0 also. This is very long-winded in python vs. IDL!

ftilt = np.array(range(len(tilt)), dtype = 'float')

for i in range(len(tilt)):
  val = tilt[i]
  if val == '': 
    val = '0'
  ftilt[i] = float(val)

is_pluto_blue = np.sum( (object == 'Pluto')     & (ftilt > 13) & (ftilt < 14) )
is_pluto_red  = np.sum( (object == 'Pluto')     & (ftilt > 20) & (ftilt < 21) )

is_hd_blue    = np.sum( (object == 'HD 146233') & (ftilt > 13) & (ftilt < 14) )
is_hd_red     = np.sum( (object == 'HD 146233') & (ftilt > 20) & (ftilt < 21) )

# Loop over every night

for jd_i in jd_uniq:

# Calc the mean sub-solar longitude for obs on this particular JD

    w = np.where( (np.trunc(jd) == jd_i) & (instrument == 'RSS') & (object == 'Pluto') & (obsmode == 'SPECTROSCOPY') )

    subsollon_mean = np.mean(subsollon[w])
    if (subsollon_mean < 0):
      subsollon_mean += 2. * math.pi

    lines.append("\n")
    lines.append('JD: ' + repr(jd_i) + ", " + date_obs[jd_i == np.trunc(jd)][0] + ", night " + repr(night_num) + 
                 ", longitude_mean_pl_rss = " + repr(subsollon_mean * r2d) )

    night_num += 1

    formatstr =(" {:<8.3f} " + # dJD 
                " {:5} " +     # UTC
		" {:>7.2f} " + # Exptime
		" {:9} " +    # Object
		" {:10} "     + # Tilt
		" {:15} " +    # Instrument / Mode
		" {:4} "     +  # Lampid
		" {:7} "     +  # grating
		" {:18} "     +  # File
		" {:>6.2f} " +  # Longitude
		" {:^4} "     +  # i
                " {:3} {:3} {:3} {:3} {:3}")

    formatstr_header = formatstr.replace("f", "")
    
    lines.append(formatstr_header.format("dJD", "UTC", "Exptime", "Object", "Tilt", "Instrument / Mode", 
                                  "Lamp", "Grating", "File", "Longit", "i", "SpI", "CR", "Rct", "Sky", "Ext"))

# Loop over each file for that night

    w = np.where((np.trunc(jd) == jd_i))[0]
    for w_i in w:
        
        have_cr       = 'Y' if len(glob.glob( files[w_i].replace('.fits', '_crmode=median_criter=3.fits'))) > 0 else ''
        have_rect     = 'Y' if len(glob.glob( files[w_i].replace('.fits', '_crmode=median_criter=3_rect.fits'))) > 0 else ''
        have_sky      = 'Y' if len(glob.glob( files[w_i].replace('.fits', '_crmode=median_criter=3_rect_y0sky=*_dysky=*.fits'))) > 0 else ''
        have_ident    = 'Y' if (call(["grep", files_short[w_i], calfile]) == 0) else ''
        have_ext      = 'Y' if len(glob.glob( files[w_i].replace('.fits', '*.txt'))) > 0 else ''
        
        line = formatstr.format(
	     jd[w_i] - amin(jd), 
	     (utc_obs[w_i])[0:5],
	     exptime[w_i], 
             object[w_i], 
	     tilt_anno[w_i], 
	     instrument[w_i] + ' / ' + obsmode[w_i], 
	     lampid[w_i], 
	     grating[w_i], 
	     files_short[w_i],
	     subsollon[w_i] * r2d, 
             w_i, 
	     have_ident, have_cr, have_rect, have_sky, have_ext)

# Finally, save the line, for the simple FITS file
    
        if ((object[w_i] not in {'FLAT'}) & 
            (proposer[w_i] == 'Throop') & 
            (files_short[w_i].find('_rect') == -1) &
            (files_short[w_i].find('_cr') == -1)):
              lines.append(line)

# Now after all the lines for this night are printed, do some statistics on the number of observations

    w_pluto_blue = where( (object == 'Pluto')     & (ftilt > 13) & (ftilt < 14))[0]
    w_pluto_red  = where( (object == 'Pluto')     & (ftilt > 20) & (ftilt < 21))[0]

    w_hd_blue    = where( (object == 'HD 146233') & (ftilt > 13) & (ftilt < 14) )[0]
    w_hd_red     = where( (object == 'HD 146233') & (ftilt > 20) & (ftilt < 21) )[0]

    w_is_raw_fits = where(is_raw_fits)[0]
    
    lines.append( "Pluto Blue: " + repr(len( set(w).intersection(w_pluto_blue).intersection(w_is_raw_fits))) )
    lines.append( "Pluto Red:  " + repr(len( set(w).intersection(w_pluto_red).intersection(w_is_raw_fits))) )

    lines.append( "HD Blue:    " + repr(len( set(w).intersection(w_hd_blue).intersection(w_is_raw_fits))) )
    lines.append( "HD Red:     " + repr(len( set(w).intersection(w_hd_red).intersection(w_is_raw_fits))) )

# Write the whole thing to a text file

np.savetxt(file_out, lines, fmt="%s")
print "Wrote: " + file_out

quit

##########
# Make a list of all of the HD and Pluto observations, along with JD
# This is so I can compute sub-solar longitude, etc... which I can't yet do in Python.
# Makes files called pluto_obs_list.txt, which I saved as txt file.
##########

for i in range(size(files)):
    if ((object[i] == 'Pluto') | (object[i] == 'HD 146233')):
        print "{0}, {1}, {2}, {3}, {4}, {5}".format(jd[i], object[i], exptime[i], tilt[i], instrument[i], obsmode[i])

# Make another useful table. This has one entry per JD. (It looks like each night spans only one JD, conveniently.)
# For each JD, it lists the UT of each Pluto spectroscopy observation.
# Then, it averages these, to give an average UT, and an average JD.
# Then I can use these to compute (in IDL) the sub-obs longitude on Pluto.
# This file is saved as 'salt_spectra_list_jd.txt'

night_num = 1

name_observer = 'Earth'
name_target   = 'Pluto'
method = 'Near Point'
fixref = 'IAU_PLUTO'
abcorr = 'LT'
r2d    = 360. / 2. / math.pi

for jd_i in jd_uniq:   # Loop over days

    print "\n"
    print 'JD: ' + repr(jd_i) + ", " + date_obs[jd_i == np.trunc(jd)][0] + ", night " + repr(night_num)
    night_num += 1

# Extract just the Pluto Spectroscopy observations
# But problem: These are actually files. We want to keep just the first .fits file, not _cr_extract.fits, etc. How??
# A: Not very elegant, but we'll just do an 'if' statement later on down.

    w = where( (obsmode == 'SPECTROSCOPY') & (object == 'Pluto') & (trunc(jd) == jd_i) ) [0]

    formatstr = " {0:<10.3f} {1:6} {2:10} {3:10} {4:12} {5:>10.2f} {6:10} {7:10} {8:10} {9:20}"

    numjd = 0
    jdsum = 0

    for w_i in range(w.size):  # Loop over observations within the day
      if ((files_short[w[w_i]].find('_rect') == -1) &
          (files_short[w[w_i]].find('_cr') == -1)):

          line = formatstr.format(jd[w[w_i]] - amin(jd), (utc_obs[w[w_i]])[0:5],
             object[w[w_i]], instrument[w[w_i]], obsmode[w[w_i]], exptime[w[w_i]], lampid[w[w_i]], grating[w[w_i]], tilt[w[w_i]], 
	     files_short[w[w_i]])

	  jdsum += jd[w[w_i]]
	  numjd += 1

#           print line

    if (numjd > 0):

# Calculate the subsolar longitude based on time
    
      et = cspice.utc2et("JD" + repr(jdsum/numjd))

# NB: Using deprecated SUBSOL method instead of SUBSLR. Mark's python library doesn't have SUBSLR.

      spoint = cspice.subsol(method, name_target, et, abcorr, name_observer)
        
      (radius, lon, lat) = cspice.reclat(spoint)

      print 'num_spect_pluto = ' + repr(numjd) + "; meanjd = " + repr(jdsum/numjd) + " for " + date_obs[jd_i == np.trunc(jd)][0] + \
          ", subsollon = " + repr(lon*r2d) + ", subsollat = " + repr(lat*r2d)

