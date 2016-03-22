# These are batch commands to process SALT data using PyRAF / PySALT
#
# Program reads all of the reduced SALT data from the SALT pipeline (.txt files).
# It them combines them into merged spectra (that is, merging red and blue tilts into one)
# And makes plots of Pluto, HD star, etc.
#
# It also takes the combined spectra and outputs them to text files.
#
# It does not do any analsis of spectral lines, etc.
#
# HBT 15-Aug-2014

# Set the prompt color. Can't figure out how to put magic command into a config file in .ipython...

# %config PromptManager.in_template = r'{color.LightGreen}In [{color.LightGreen}{count}{color.LightGreen}]: '

#from pyraf import iraf
#from iraf import pysalt
#from iraf import saltspec
#from iraf import specextract
#from iraf import specprepare
#from iraf import specidentify
#from iraf import specrectify
#from iraf import specsky

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

def gaus(x,a,x0,y0, sigma):
  return y0 + a*exp(-(x-x0)**2/(2*sigma**2))

##########
# Remove artifacts from the chip gaps
##########
#
# Chip edges. I could design a fancy algorithm. But looking at the data, it looks like the gaps
# are like this:
#   [GOOD DATA] [2-4 PIXELS BAD DATA] [~100 0's IN GAP]  [2-4 PIXELS BAD DATA] [GOOD DATA]
#
# Same at the edges too. Remove 4 pixels.
#
# In some cases the bad data is only 2 pixels wide, or at the far edges even less.
# But in the general case it's 4, so we remove all of those.
#
# I'd like this routine to do CR subraction too. But that is harder.

def clean_artifacts(arr):

  width_gap_edges = 4
  is_gap = (arr == 0)		# Get the gap portions

# Convolve with a boxcar kernel. This is to widen the gaps. Convolving with a 9-element kernel makes 
# array 4 elements longer on either end... remove this.

  is_gap_wider = (int_(convolve(is_gap, Box1DKernel(width_gap_edges*2+1))==0))[width_gap_edges:-width_gap_edges]  == 0

  is_ends = np.zeros(is_gap.size)
  is_ends[0:width_gap_edges] = True
  is_ends[-width_gap_edges:] = True

  is_bad = (np.array(is_gap_wider) + np.array(is_ends)) > 0
  is_good = (is_bad == 0)

  arr_clean = np.copy(arr)
  arr_clean[is_bad] = 0
  return arr_clean

##########
# Function to merge multiple (>=1) spectra from red and/or blue into one
##########

def merge_multiple_rss(fluxes_blue_in, wavelengths_blue_in, 
                       fluxes_red_in, wavelengths_red_in, 
                       wavelength_merged):
#   
#                          flux_blue_in_0, wavelength_blue_in_0,
#                          flux_red_in_0, wavelength_red_in_0,
#                          flux_blue_in_1, wavelength_blue_in_1,
#                          flux_red_in_1, wavelength_red_in_1,

# First resample each spectrum into a new common baseline -- that is, a standard evenly gridded wavelength

    num_red_in  = shape(fluxes_red_in)[0]		# Number of red spectra passed
    num_blue_in = shape(fluxes_blue_in)[0]		# Number of blue spectra passed
    length_spectrum = shape(wavelength_merged)[0]		# Shape of output drizzled spectrum

# Create the output where the drizzled fluxes go

    fluxes_red  = np.zeros((num_red_in,  length_spectrum), dtype = 'f')
    fluxes_blue = np.zeros((num_blue_in, length_spectrum), dtype = 'f')

    for i in range(num_red_in):
      print 'Drizzling red ' + repr(i)
      fluxes_red[i,:] = drizzle_flux_wavelength(fluxes_red_in[i,:], wavelengths_red_in[i,:], wavelength_merged)

    for i in range(num_blue_in):
      print 'Drizzling blue ' + repr(i)
      fluxes_blue[i,:] = drizzle_flux_wavelength(fluxes_blue_in[i,:], wavelengths_blue_in[i,:], wavelength_merged)
   
#       flux_red_0  = drizzle_flux_wavelength(flux_red_in_0,  wavelength_red_in_0,  wavelength_merged)
#       flux_blue_1 = drizzle_flux_wavelength(flux_blue_in_1, wavelength_blue_in_1, wavelength_merged)
#       flux_red_1  = drizzle_flux_wavelength(flux_red_in_1,  wavelength_red_in_1,  wavelength_merged)

# Now really merge the four spectra into one. 
# This does require rescaling the spectra... because of SALT's variable pupil size, spectra at 
# 'identical' settings will have dfft amplitudes
#
# The funky int_ is a divisor: it is equal to 1 if only one of (red0, red1) is set... otherwise, if both are set, it's set to 2, so we divide by 2.
# This is the right way to combine two signals: take average if they're both there, or just one if only one is there.

# Now combine all the red spectra
# *** Ahh, we can't just combine these as an arithmetic mean. We have to weight them, because of SALT's variable pupil size.

    num_vals_red  = np.zeros(length_spectrum)
    num_vals_blue = np.zeros(length_spectrum)

    flux_red  = np.zeros(length_spectrum)
    flux_blue = np.zeros(length_spectrum)

#       print "num_vals_red =" + repr(num_vals_red)
#       print "num_vals_blue =" + repr(num_vals_blue)

# At each wavelength point, count how many valid red flux data point we have

    for i in range(num_red_in):
      num_vals_red += (fluxes_red[i,:] > 0).astype(int)

    indices_common_red = where(num_vals_red == num_red_in)		# Get indices of the wavelengths which are common between all spectra
    weight_red = np.zeros(num_red_in)
    weight_blue = np.zeros(num_red_in)

# For each spectrum, gets the median value of the flux at the shared wavelengths, and save it

    for i in range(num_red_in):
      weight_red[i] = median(fluxes_red[i,indices_common_red])
    weight_red /= np.max(weight_red)					# Normalize it

# Count how many valid blue flux data points we have

    for i in range(num_blue_in):
      num_vals_blue += (fluxes_blue[i,:] > 0).astype(int)

    indices_common_blue = where(num_vals_blue == num_blue_in)
    for i in range(num_blue_in):
      weight_blue[i] = median(fluxes_blue[i,indices_common_blue])
    weight_blue /= np.max(weight_blue)					# Normalize it

# Sum all the individual spectra. Weight each wavelength by the median of that spectrum, and by the number of spectra that observe that wavelength

    for i in range(num_red_in):
      flux_red += (fluxes_red[i,:]) / weight_red[i]
    flux_red /= num_vals_red

    for i in range(num_blue_in):
      flux_blue += (fluxes_blue[i,:]) / weight_blue[i]
    flux_blue /= num_vals_blue

# Get rid of the NaNs (which are caused by 0/0 division, for wavelengths that have no flux)

    flux_red[num_vals_red==0] = 0
    flux_blue[num_vals_blue==0] = 0

# Compute the ratio blue:red in the overlap region

    w_overlap = (flux_blue * flux_red) > 0 	# Get indices where blue >0 and red>0 (ie, the overlap portion)
    ratio_b_r = np.median( (flux_blue / flux_red)[w_overlap] )

#       print "w_overlap = " + repr(w_overlap)
#     
#       print "ratio_b_r = " + repr(ratio_b_r)

# Actually do the merging itself

    flux_merged = flux_blue + ratio_b_r * flux_red

    flux_merged[w_overlap] /= 2

# Make some diagnostic plots to confirm
   
#       plot(wavelength_merged, 1.0 * flux_red / np.median(flux_red[w_overlap]))
#       plot(wavelength_merged, 1.2 * flux_blue / np.median(flux_red[w_overlap]))
#       plot(wavelength_merged, 1.5 * flux_merged / np.median(flux_merged[w_overlap]))
#       plt.title('Merged')
#       show()

# Normalize the data by its median value

    flux_merged = flux_merged * 0.5 / np.median(flux_merged)

    return flux_merged
  
##########
# Function to merge four spectra (2x red, 2x blue) into one
##########

def merge_four_rss(flux_blue_in_0, wavelength_blue_in_0,
                flux_red_in_0, wavelength_red_in_0,
                flux_blue_in_1, wavelength_blue_in_1,
                flux_red_in_1, wavelength_red_in_1,
                wavelength_merged):

# First resample each spectrum into a new common baseline -- that is, a standard evenly gridded wavelength

    flux_blue_0 = drizzle_flux_wavelength(flux_blue_in_0, wavelength_blue_in_0, wavelength_merged)
    flux_red_0  = drizzle_flux_wavelength(flux_red_in_0,  wavelength_red_in_0,  wavelength_merged)
    flux_blue_1 = drizzle_flux_wavelength(flux_blue_in_1, wavelength_blue_in_1, wavelength_merged)
    flux_red_1  = drizzle_flux_wavelength(flux_red_in_1,  wavelength_red_in_1,  wavelength_merged)

# Now really merge the four spectra into one. 
# This does require rescaling the spectra... because of SALT's variable pupil size, spectra at 
# 'identical' settings will have dfft amplitudes
#
# The funky int_ is a divisor: it is equal to 1 if only one of (red0, red1) is set... otherwise, if both are set, it's set to 2, so we divide by 2.
# This is the right way to combine two signals: take average if they're both there, or just one if only one is there.
    
    flux_blue = (flux_blue_0 + flux_blue_1) / (1 + int_(flux_blue_0 * flux_blue_1 > 0))
    flux_red  = (flux_red_0  + flux_red_1)  / (1 + int_(flux_red_0  * flux_red_1  > 0))

# Compute the ratio blue:red in the overlap region

    w_overlap = (flux_blue * flux_red) > 0 	# Get indices where blue >0 and red>0 (ie, the overlap portion)
    ratio_b_r = np.median( (flux_blue / flux_red)[w_overlap] )

#     print "ratio_b_r = " + repr(ratio_b_r)

# Actually do the merging itself

    flux_merged = flux_blue + ratio_b_r * flux_red
    flux_merged[w_overlap] /= 2

# Normalize the data by its median value

    flux_merged = flux_merged * 0.5 / np.median(flux_merged)
 
#     plot(flux_merged)
#     plt.title('flux_merged')
#     show()
# 
#     plot(wavelength_merged, flux_blue_0)
#     plt.title('flux_blue_0')
#     show()
# 
#     plot(wavelength_merged, flux_blue_1)
#     plt.title('flux_blue_1')
#     show()
# 
#     plot(wavelength_merged, flux_red_0)
#     plt.title('flux_red_0')
#     show()
# 
#     plot(wavelength_merged, flux_red_1)
#     plt.title('flux_red_1')
#     show()
# 
#     plot(wavelength_merged, flux_red)
#     plt.title('flux_red')
#     show()
# 
#     plot(wavelength_merged, w_overlap)
#     plt.title('w_overlap')
#     plt.ylim((-1,2))
#     show()
# 
#     plot(wavelength_merged, flux_merged)
#     plt.title('w_overlap')
#     plt.ylim((-1,2))
#     show()

    return flux_merged

##########
# Resample a spectrum onto a standard evenly gridded wavelength array    
##########
# *** I am sure my code is very slow, and this routine is a real bottleneck!
##########
# *** Also, this routine really badly kills SNR. If I take an 0.8A input, and resample to 0.5A output, 
#     it just copies one value (the closest bin) -- rather than resampling them all. I should essentially be doing
#     a boxcar smoothing on this data *before* I resample it, at least if I resample to lower res.
#     (If I resample to 5x nominal, I should use a boxcar smoothing width of 5, I think.)
#     
#     NB: Drizzle changes both the scale and the offset, so I need to do the same here too.
#    


def drizzle_flux_wavelength(flux_in, wavelength_in, wavelength_out):

    dw_native = wavelength_in[1] - wavelength_in[0]   # Current resolution
    dw_out    = wavelength_out[1] - wavelength_out[0] # Output resolution

    flux_out = np.empty_like(wavelength_out)

    DO_METHOD_ZOOM = False

    if (DO_METHOD_ZOOM):

# Zoom it

# Now roll it so that it lines up with wavelength_out. Problem is that wavelength_in and wavelength_out have different ranges.
# I could just build a new array 

        range_alam_in  = (np.amin(wavelength_in),  np.amax(wavelength_in))
        range_alam_out = (np.amin(wavelength_out), np.amax(wavelength_out))

        delta_alam_in  = wavelength_in[1]  - wavelength_in[0]
        delta_alam_out = wavelength_out[1] - wavelength_out[0]

# We want to take the input array, and pad it on both ends to make it match the limits of the output array.
# Then when we apply the .zoom(), it will all scale properly to the output wavelength array.

# Create the padded version of the input array

        sizex_flux_in_padded = np.round( (range_alam_out[1] - range_alam_out[0]) / delta_alam_in) 

        print "flux_in.shape = " + repr(flux_in.shape)
        flux_in_padded = np.zeros(sizex_flux_in_padded)
        print "flux_in_padded.shape = " + repr(flux_in_padded.shape)
        index_start = np.round( (wavelength_in[0] - wavelength_out[0]) / (delta_alam_in) ) 
        flux_in_padded[index_start:index_start + len(flux_in)] = flux_in

        zoom = len(flux_out) * 1. / len(flux_in_padded)
        print "zoom2 = " + repr(zoom)

        flux_out_test = scipy.ndimage.interpolation.zoom(flux_in_padded, zoom)
        print "flux_out_test.shape = " + repr(flux_out_test.shape)
        print "flux_out_test.shape = " + repr(flux_out_test.shape)
        print "wavelength_out.shape = " + repr(wavelength_out.shape)

        return flux_out_test

    bin_max = flux_in.size - 1 # Index of the highest bin in the input array
  
    # Apply a boxcar smoothing to it
    # Philosophy here is that we are downsampling here to a common resolution. 
    # We only want to apply this here if we are really downsampling. In most cases we
    # don't want to do this -- we want to preserve all resolution into the output text files,
    # and then downsample there.

#       flux_in_smoothed = convolve(flux_in, Box1DKernel(dw_out / dw_native))

    for i in range(wavelength_out.size):
        d_w = abs(wavelength_in - wavelength_out[i])
        bin_closest = (where(d_w == amin(d_w)))[0][0] # Both zeroes are necessary here: first is to unwrap results of where(), 
                                                      # and second is in case two bins are equal

						# If we are requesting a wavelength beyond the range of this input spectrum, set output value to 0

#           flux_out[i] = 0 if ((bin_closest == 0) | (bin_closest == bin_max)) else flux_in_smoothed[bin_closest]
        flux_out[i] = 0 if ((bin_closest == 0) | (bin_closest == bin_max)) else flux_in[bin_closest]

#       print "wavelength_in:  range " + repr(np.amin(wavelength_in)) + " .. " + repr(np.amax(wavelength_in)) + " AA, " + repr(len(wavelength_in)) + " bins"
#       print "wavelength_out: range " + repr(np.amin(wavelength_out)) + " .. " + repr(np.amax(wavelength_out)) + " AA, " + repr(len(wavelength_in)) + " bins"
#   
#       plot(wavelength_in, flux_in, color = 'blue')
#       plot(wavelength_out, flux_out + 0.1, color='red')
#       plot(wavelength_out, flux_out_test + 0.2, ls = '-', color = 'green')
#       plt.show()
#   
#       plot(wavelength_out, flux_out_test + 0.2, ls = '-', color = 'green')
#       plt.show()

    return flux_out

##########
# Main Program
##########

colors = ['red', 'blue', 'green', 'purple', 'orange', 'black', 'darkblue', 'salmon', 'pink', 'brown', 'olive', 'violet', 'brown', 'tan']
    
dir           = '/Users/throop/python'

# CD into python. Just dump everything there for now.
# I don't know how PyRAF handles directories, so better to cd into the right one from the start.

os.chdir(dir)

# Startup SPICE. We use it for time conversions only

file_kernels = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'
cspice.furnsh(file_kernels)

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

# Get a list of raw unprocessed FITS files. Do not include *_rect_sky.fits, etc.

# file_list = glob.glob(dir_data + '/mbxgpP2014*_crmode=fast_criter=3_rect_y0sky*dysky=20_dyspect=20_y0spect*.txt')

# Define a combined wavelength grid ('wavelength_merged') for all of the output spectra (HD + Pluto)

dw = 0.5    # Width of output wavelength bins in Angstrom. Original source files are 0.94 A I think.
            # In general we don't want to reduce resolution in the current algorithm: retain everything we can.

year_obs = 2015

if (year_obs == 2014):
  dir_data      = "/Users/throop/Data/SALT_Pluto_2014/product"
  file_list = glob.glob(dir_data + '/mbxgpP2014*[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits')

if (year_obs == 2015):
  dir_data      = "/Users/throop/Data/SALT_Pluto_2015/product"
  file_list = glob.glob(dir_data + '/mbxgpP2015*[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits')

# Define the parameters for the processing used on this run. This can be changed easily. All this stuff replaces the ".fits" string

extensions_processing = '_crmode=median_criter=3_rect_y0sky*dysky=20_dyspect=20_y0spect*.txt'

files = np.array(file_list)

# Read the JD from each file. Then sort the files based on JD.

jd = []

print "Reading all FITS headers for JD: " + repr(files.size) + " files"

for file in files:
    sys.stdout.write('.') # print w/o newline
    hdulist = fits.open(file)
    jd.append(hdulist[0].header['JD'])
    hdulist.close()
print

print "Sorting by JD\n"
    
indices = np.argsort(jd)
files = files[indices] # Resort by JD. ** This takes all the flats (which are alphabetically at end) and distributes them throughout
# pdb.set_trace()
print

# Now read in all the files a second time, this time in order by JD
print "Reading all FITS headers: " + repr(files.size) + " files"
  
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

print

# Convert from lists to NumPy arrays as needed

jd          = np.array(fits_jd)
utc         = np.array(fits_utc_obs)
date        = np.array(fits_date_obs)
date_evening= np.array(fits_date_obs)
airmass     = np.array(fits_airmass)
object      = np.array(fits_object, dtype='S30') # np.array can use string arrays as easily as float arrays
instrument  = np.array(fits_instrume)
proposer    = np.array(fits_proposer)
exptime     = np.array(fits_exptime)
lampid	    = np.array(fits_lampid)
grating     = np.array(fits_grating)
tilt        = np.array(fits_gr_angle)
obsmode     = np.array(fits_obsmode)  # Remove extra spaces
# pdb.set_trace()

# 2014 uses 'HD 146233', while 2015 is 'HD146233'. Probably my error in making Step-2? Make them match.

for i in range(size(object)):
  oi = object[i]
  if oi == "HD146233": object[i] = "HD 146233"

# Convert the dates. The deal is that we want to tag these observations not by the 
# actual SAST calendar date, but by the date of the evening that the observation started.

for i in range(len(date_evening)):
    s    = cspice.et2utc(cspice.utc2et('JD' + repr(trunc(jd[i]))), 'C', 3) # JD starts at noon UT! So we trunc() the JD.
    strs = s.split(" ")						  # SPICE returns 'OCT' and we must convert to 08.
    									  # SPICE will not do conversion itself.
    month_num = { 'JAN' : '01', 'FEB' : '02', 'MAR' : '03', 'APR' : '04', 'MAY' : '05', 'JUN' : '06', 
                  'JUL' : '07', 'AUG' : '08', 'SEP' : '09', 'OCT' : '10', 'NOV' : '11', 'DEC' : '12' }[strs[1]]
    date_evening[i] = strs[0] + '-' + month_num + '-' + strs[2]

# Extract just the good ones
# ** When initializing the np string array, need to set the maximum string length.
# ** I choose a long default here, so we won't exceed it.

# files      = np.array(file_list,dtype='S100') # Was bug here: I was copying file_list (not JD_sorted), not files (is JD sorted)

files_short = np.array(files)

for i in range(files.size):
    files_short[i] = files[i].split('/')[-1]  # Get just the filename itself

# Filter out the observations which are obviously not relevant... like, SALTICAM. Keep the ARCs.
   
is_bad =  ((proposer != 'Throop') | (instrument != 'RSS') | (object == 'BIAS') | (obsmode == 'FABRY-PEROT') | (obsmode == 'IMAGING') | 

   # List all bad 2014 files

          (files_short == "mbxgpP201408260016.fits") | # Bad ARC... all noise. Don't know why. It was retaken, so skip it.
          (files_short == "mbxgpP201408160094.fits") | # Bad HD red. Took for 2 sec immediately after good 4-sec with same settings. 
                                                       # Took 3 not 2.
          (files_short == "mbxgpP201409050029.fits") | # Took four HD Red instead of two. Don't know why. Deleting first two.
          (files_short == "mbxgpP201409050030.fits") | # Took four HD Red instead of two. Don't know why. Deleting first two.
          (files_short == "mbxgpP201409050198.fits") | # Took four Pluto red instead of two. Don't know why. Deleting first two
          (files_short == "mbxgpP201409050199.fits") | # Took four Pluto red instead of two. Don't know why. Deleting first two
          (files_short == "mbxgpP201409050189.fits") |   # Took 3 Pluto blue not 2. From order of obs, maybe this was off-target? 

    # List all bad 2015 files

          (files_short == "mbxgpP201505040066.fits") |   # Took 3 Pluto red not 2. This one is OK, but off target and not extracted. ** Should use this one!
#             (files_short == "mbxgpP201505040064.fits") |   # Bad interference fringes. Can use, but poor beyond 7600
#             (files_short == "mbxgpP201505040065.fits") |   # Bad interference fringes. Can use, but poor beyond 7600
	  (files_short == "mbxgpP201506290031_ORIGINAL") | # SALT error in keywords -- I fixed but keep original
	  (files_short == "mbxgpP201506290032_ORIGINAL") | # SALT error in keywords -- I fixed but keep original
	  (files_short == "mbxgpP201506290028.fits") | # Operator error -- this is proper 1-sec HD; was followed w 30-sec HD which I use
	  (files_short == "mbxgpP201506290036") |     # Operator error? Took 7 HD Red. All look OK. So I keep the 4 earliest and
	  (files_short == "mbxgpP201506290038") |     # longest ones (because arcs were not done immediately after all tilts)
	  (files_short == "mbxgpP201506290039") |     # and tossed remaining 3. Should read nightlog to get full story.
	  (files_short == "mbxgpP201506300027.fits") |     # SA took five. First four are OK. Last one is off chip. Could fix, but toss it.
	  (files_short == "mbxgpP201505030056.fits") |     # Bad SPECRECTIFY (from arc in 0053). Wrote to salthelp 1-Nov-2015
	  (files_short == "mbxgpP201505030057.fits") |     # Bad SPECRECTIFY (from arc in 0053). Wrote to salthelp 1-Nov-2015
	  (files_short == "mbxgpP201506300001.fits"))      # _RECT looks weird. SA took two. I keep the first and toss second.

          

# is_bad = (proposer != 'Throop') | (instrument != 'RSS') | (object == 'BIAS') | (obsmode == 'FABRY-PEROT') | (obsmode == 'IMAGING')
        
is_good = (is_bad == False)

print "Removing other PIs, Salticam, F-P, BIAS: down-selecting " + repr(files.size) + " to " + repr(sum(is_good))

files_copy = files.copy()
files      = files_copy[is_good]

# files      = files[is_good]  # Is this legal in python? Does it make sense, or is it a view into itself??
files_short=files_short[is_good]
jd         = jd[is_good]
utc        = utc[is_good]
date       = date[is_good]
date_evening = date_evening[is_good]
airmass    = airmass[is_good]
object     = object[is_good]
instrument = instrument[is_good]
proposer   = proposer[is_good]
exptime    = exptime[is_good]
lampid     = lampid[is_good]
grating    = grating[is_good]
obsmode    = obsmode[is_good]
tilt       = tilt[is_good].astype(float) # For some reason this reads as a string from the FITS file...

# Make a NP-compatible copy of 'files' as well... as array, not list!

# pdb.set_trace()

files_rect         = np.zeros(size(files),dtype='S200')	# Originally I was using string.copy. But that can change dtype to shorter!
files_rect_sky     = np.zeros(size(files),dtype='S200')
files_rect_sky_ext = np.zeros(size(files),dtype='S200')

files_cr_rect         = np.zeros(size(files),dtype='S200')
files_cr_rect_sky     = np.zeros(size(files),dtype='S200')
files_cr_rect_sky_ext = np.zeros(size(files),dtype='S200')

for i in range(files.size):
    files_rect[i]         = string.replace(files[i], '.fits', '_rect.fits')
    files_rect_sky[i]     = string.replace(files[i], '.fits', '_rect_sky.fits') 
    files_rect_sky_ext[i] = string.replace(files[i], '.fits', '_rect_sky_ext.txt')

    files_cr_rect[i]         = string.replace(files[i], '.fits', '_cr_rect.fits')
    files_cr_rect_sky[i]     = string.replace(files[i], '.fits', '_cr_rect_sky.fits') 
    files_cr_rect_sky_ext[i] = string.replace(files[i], '.fits', '_cr_rect_sky_ext.txt')

# Make arrays for the various output files

# Get a list of Arcs, Flats, etc. Every image is exactly one of these.

print "Filtering by object name\n"

is_arc   = (object == 'ARC')
is_pluto = (object == 'Pluto')
is_hd    = (object == 'HD 146233')
is_flat  = (object == 'FLAT')

calfile = "salt_pluto_cal.txt"

##########
# Now plot the spectra
##########

# First find an appropriate spectrum (.txt) file for each .fits file
# If a file is missing, that is probably OK -- e.g., it might be an ARC, which makes no .txt file

files_spec = np.zeros(size(files), dtype='S150')

print "Matching .fits file to .txt files"

for i in range(size(files_spec)):
    s = str.replace(files[i], '.fits', extensions_processing)
    g = glob.glob(s)			# Find a .txt file that matches
    if (size(g) > 0):
        g0 = g[0].split('/')[-1]		# Get just the filename itself
        print files_short[i] + " -> " + g0
        files_cr_rect_sky_ext[i] = g0
    else:
        print files_short[i] + " : No .txt file found. Object = " + object[i]

 # pdb.set_trace()

fluxes        = []
dfluxes       = []
wavelengths   = []

files_spec = files_cr_rect_sky_ext

# *** Should look up 'genfromtext' or 'getfromcsv' in numpy
# Actually, SciPy FAQ says that numpy.loadtxt() is probably the way to do it.
# It will deal with header lines, CSV, etc.
# Astropy also has a 'unified table reader'

# Read all of the data into a list of arrays.
# It might be possible to do this as an array of arrays -- not sure.
# See 'advanced topics' section in NumPy book.
# File is a simple three-column list of numbers, separated by spaces

#for file in files_spec: 
#    d = loadtxt(file) # Oh my god! This works like magic.
#    wavelengths.append(d[:,0])
#    fluxes.append(d[:,1])
#    dfluxes.append(d[:,2])

rcParams['figure.figsize'] = 15, 10 # Make plot normal size

#colors = np.array(['grey', 'brown', 'orange', 'yellow', 'blue', 'red', 'pink', 'green'])
#for i in range(len(fluxes)):
#    f = fluxes[i].astype(float)
#    f = f / np.amax(f)  # normalize it
#    plot(wavelengths[i],f+i, color=colors[i % colors.size], ls='-')
#    plt.axis([3000,9000,0,100])
  
# plt.xlim((3000,10000))
    
##########
# Now get the indices for the exposures for the HD spectra
##########

w_spect_hd_blue  = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 13) & (tilt < 14))[0] # length 10, for 2014
w_spect_hd_red   = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 20) & (tilt < 21))[0] # Length 10, for 2014

##########
# Now get the indices for the exposures for the Pluto spectra
##########

w_spect_pluto_blue  = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 13) & (tilt < 14))[0] # Length 27, for 2014
w_spect_pluto_red = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 20) & (tilt < 21))[0]   # Length 28, for 2014

wavelengths_hd_blue = []  # 'wavelengths', 'fluxes', 'dfluxes' are lists, not np arrays, because they are flexible length.
fluxes_hd_blue = []
dfluxes_hd_blue = []
wavelengths_hd_red  = []
fluxes_hd_red  = []
dfluxes_hd_red  = []

wavelengths_pluto_blue = []
fluxes_pluto_blue = []
dfluxes_pluto_blue = []
wavelengths_pluto_red  = []
fluxes_pluto_red  = []
dfluxes_pluto_red  = []

# Read in the HD spectra from disk. Up til now, everything has been using FITS headers.

rcParams['figure.figsize'] = 15,10 # Make plot normal size

#quit

for i in w_spect_hd_red: # All good
    file = files_cr_rect_sky_ext[i]
    d = loadtxt(dir_data + '/' + file)
    d[d<0] = 0		# Clamp the values. A few of the txt files from SPECEXTRACT have negatives around the chip gaps -- apparently a bug
    wavelengths_hd_red.append(d[:,0])
    fluxes_hd_red.append(d[:,1])
    dfluxes_hd_red.append(d[:,2]) 

for i in w_spect_hd_blue:  # one missing file. No 0032 file. 
    file = files_cr_rect_sky_ext[i]
    d = loadtxt(dir_data + '/' + file)
    d[d<0] = 0
    wavelengths_hd_blue.append(d[:,0])
    fluxes_hd_blue.append(d[:,1])
    dfluxes_hd_blue.append(d[:,2])     

# Read in the Pluto spectra

for i in w_spect_pluto_red: # All good
    file = files_cr_rect_sky_ext[i]
    d = loadtxt(dir_data + '/' + file)
    d[d<0] = 0
    wavelengths_pluto_red.append(d[:,0])
    fluxes_pluto_red.append(d[:,1])
    dfluxes_pluto_red.append(d[:,2]) 

for i in w_spect_pluto_blue:  # One missing file
    file = files_cr_rect_sky_ext[i]
    d = loadtxt(dir_data + '/' + file)
    d[d<0] = 0
    print "Reading Pluto Blue file " + file
    wavelengths_pluto_blue.append(d[:,0])
    fluxes_pluto_blue.append(d[:,1])
    dfluxes_pluto_blue.append(d[:,2])    
    
# Set the plot size
    
rcParams['figure.figsize'] = 15, 10

######
# Make a plot of Ha and OII before the shifting
######

DO_SHIFT = False

range_wavelength = 100
wavelength_band = np.array((6563, 7605))  # O II
name_band = np.array(("Ha", "OII"))
ylim_max = np.array((25000,20000))

if (DO_SHIFT):

  for j in range(name_band.size):

    for i in range(w_spect_pluto_red.size):
        plot(wavelengths_pluto_red[i], fluxes_pluto_red[i] + i*0.1, color=colors[i % len(colors) ])
        bin = wheremin(abs(wavelengths_pluto_red[i] - wavelength_band[j]))
        plt.text(wavelengths_pluto_red[i][bin] + 30, fluxes_pluto_red[i][bin], 
	files_short[w_spect_pluto_red[i]], color=colors[i % len(colors)] ) 
   
    plt.title('Pluto @ ' + name_band[j] + ' =' + repr(wavelength_band[j]) + " $\AA$ [not shifted]", fontsize=24)
    plt.xlim((-range_wavelength/2 + wavelength_band[j], range_wavelength/2 + wavelength_band[j]))
    plt.ylabel('Intensity [arbitrary]', fontsize=24)
    plt.xlabel('Wavelength [$\AA$]', fontsize=24)
    plt.ylim((0, ylim_max[j]))

    show()

#stop

# Do some wavelength calibration. 
# Take all of the red-tilt spectra, and shift them to make sure they match H-alpha (6562.8 angstroms). 
# *** For some reason this doesn't seem to be affecting the output files properly. I'm not sure why.
#     Also, there may be a methodology issue: even if Ha is shifted to match, then O2 @ 7620 doesn't match. 
#     So maybe it's hopeless.
#     Also, this really points to a problem with how SPECIDENTIFY has been run. Better to fix it there, than here.

# Shift based on Ha. Then plot both Ha and OII again.

  index_extract = np.array((2,2))
  wavelength_ha = wavelength_band[0]

#stop

  wavelength_reference = wavelength_band[1]  # Which should we use as a reference wavelength: Ha, OII, etc.

  j=1
  p0_oii = [-30000,7605,35000,5]

  wavelength_shift_oii = np.zeros(w_spect_pluto_red.size)

  for i in range(w_spect_pluto_red.size):
    data_x = wavelengths_pluto_red[i][1510:1540]
    data_y = fluxes_pluto_red[0][1510:1540]
    popt,pcov = curve_fit(gaus, data_x, data_y, p0=p0_oii)
    wavelength_shift_oii[i] = 7605 - popt[1]
    print "Center at " + repr(popt[1]) + "; shift = " + repr(wavelength_shift_oii[i])

  bin_shift_oii =  -int_(np.round(wavelength_shift_oii / (wavelengths_pluto_red[1][1]-wavelengths_pluto_red[0][0])))

# Make a one-off plot of OII band, shifted

  for i in range(w_spect_pluto_red.size):
    plot(wavelengths_pluto_red[i], roll(fluxes_pluto_red[i] + i*0.1, bin_shift_oii[i]), color=colors[i % len(colors) ])
    bin = wheremin(abs(wavelengths_pluto_red[i] - wavelength_band[1]))
    plt.text(wavelengths_pluto_red[i][bin] + 30, fluxes_pluto_red[i][bin], 
    files_short[w_spect_pluto_red[i]], color=colors[i % len(colors)] ) 

  plt.title('Pluto @ ' + name_band[j] + ' =' + repr(wavelength_band[j]) + " $\AA$ [shifted]", fontsize=24)
  plt.xlim((-range_wavelength/2 + wavelength_band[j], range_wavelength/2 + wavelength_band[j]))
  plt.ylabel('Intensity [arbitrary]', fontsize=24)
  plt.xlabel('Wavelength [$\AA$]', fontsize=24)
  plt.ylim((0, ylim_max[j]))
  
  show()

# Calculate the shift for Ha @ 6563 AA

  p0_ha = [-10000,6563,25000,5]
  j=0

  wavelength_shift_ha = np.zeros(w_spect_pluto_red.size)

  for i in range(w_spect_pluto_red.size):
    data_x = wavelengths_pluto_red[i][350:460]
    data_y = fluxes_pluto_red[0][350:460]
    popt,pcov = curve_fit(gaus, data_x, data_y, p0=p0_ha)
    wavelength_shift_ha[i] = 6563 - popt[1]
    print "Center at " + repr(popt[1]) + "; shift = " + repr(wavelength_shift_ha[i])

  bin_shift_ha =  int_(np.round(wavelength_shift_ha / (wavelengths_pluto_red[1][1]-wavelengths_pluto_red[0][0])))

# NB: This shift has a different sign than the one above! In calculating bin_shift_ha

# Make a one-off plot of Ha band, shifted

  for i in range(w_spect_pluto_red.size):
    plot(wavelengths_pluto_red[i], roll(fluxes_pluto_red[i] + i*0.1, bin_shift_ha[i]), color=colors[i % len(colors) ])
    bin = wheremin(abs(wavelengths_pluto_red[i] - wavelength_band[1]))
    plt.text(wavelengths_pluto_red[i][bin] + 30, fluxes_pluto_red[i][bin], 
    files_short[w_spect_pluto_red[i]], color=colors[i % len(colors)] ) 

  plt.title('Pluto @ ' + name_band[j] + ' =' + repr(wavelength_band[j]) + " $\AA$ [shifted]", fontsize=24)
  plt.xlim((-range_wavelength/2 + wavelength_band[j], range_wavelength/2 + wavelength_band[j]))
  plt.ylabel('Intensity [arbitrary]', fontsize=24)
  plt.xlabel('Wavelength [$\AA$]', fontsize=24)
  plt.ylim((0, ylim_max[j]))

  show()

colors = ['red', 'blue', 'green', 'purple', 'orange', 'black', 
          'darkblue', 'salmon', 'pink', 'brown', 'olive', 'violet', 'brown', 'tan']

##########
# Make a plot of HD
# For each date, show all of the data observations (not arcs)
##########

dy = 3.5						# Separation between y lines
y_0 = 0							# Current y position for this line

for i in range(w_spect_hd_blue.size):    

    if (date_evening[w_spect_hd_blue[i]] != date_evening[w_spect_hd_blue[i-1]]):	# If we're on a new date, jump to the next vertical position
        y_0 += dy
    
    plot(wavelengths_hd_blue[i], fluxes_hd_blue[i] / amax(fluxes_hd_blue[i]) + y_0, color=colors[i % len(colors)])

y_0 = 0
for i in range(w_spect_hd_red.size):    
    if (date_evening[w_spect_hd_red[i]] != date_evening[w_spect_hd_red[i-1]]):
        y_0 += dy
    plot(wavelengths_hd_red[i], fluxes_hd_red[i]   / amax(fluxes_hd_red[i])  + y_0, color=colors[i % len(colors)])
    plt.text(9400, y_0, date_evening[w_spect_hd_red[i]])

plt.ylim((0,w_spect_hd_red.size))
plt.xlabel('Wavelength [$\mu$ m]', fontsize = 20)
plt.ylabel('I [arbitrary]', fontsize = 20)
plt.title('HD', fontsize=20)
show()

#   for i in range(w_spect_hd_blue.size):    
#       plot(wavelengths_hd_blue[i], fluxes_hd_blue[i] / amax(fluxes_hd_blue[i]) + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
#       plot(wavelengths_hd_red[i], fluxes_hd_red[i]   / amax(fluxes_hd_red[i])  + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
#       plt.text(9400, i-i%2 + 0.3, date_evening[w_spect_hd_blue[i]])
#   plt.ylim((0,(w_spect_hd_red.size)/2-1))
#   plt.xlim((3000,10000))
#   plt.xlabel('Wavelength [$\AA$]', fontsize = 20)
#   plt.ylabel('I [arbitrary]', fontsize = 20)
#   plt.title('HD', fontsize=20)
#   show()

##########
# Make a plot of Pluto
# For each date, show all of the data observations (not arcs)
##########

dy = 1.4						# Separation between y lines
y_0 = 0	- 0.9 * dy						# Current y position for this line

for i in range(w_spect_pluto_blue.size):    

    if (date_evening[w_spect_pluto_blue[i]] != date_evening[w_spect_pluto_blue[i-1]]):	# If we're on a new date, jump to the next vertical position
        y_0 += dy
    
    plot(wavelengths_pluto_blue[i], fluxes_pluto_blue[i] / amax(fluxes_pluto_blue[i]) + y_0, color=colors[i % len(colors)])

y_0 = 0 - 0.9 * dy

for i in range(w_spect_pluto_red.size):    
    if (date_evening[w_spect_pluto_red[i]] != date_evening[w_spect_pluto_red[i-1]]):
        y_0 += dy
    plot(wavelengths_pluto_red[i], fluxes_pluto_red[i]   / amax(fluxes_pluto_red[i])  + y_0, color=colors[i % len(colors)])
    plt.text(9400, y_0, date_evening[w_spect_pluto_red[i]])

plt.ylim((0,11))
plt.xlabel('Wavelength [$\mu$ m]', fontsize = 20)
plt.ylabel('I [arbitrary]', fontsize = 20)
plt.title('Pluto', fontsize=20)
show()

##########
# Now clean the spectra, taking out any artifacts near chip gap edges, etc.
##########


fluxes_pluto_red_clean  = np.copy(fluxes_pluto_red)
fluxes_pluto_blue_clean = np.copy(fluxes_pluto_blue)
fluxes_hd_red_clean     = np.copy(fluxes_hd_red)
fluxes_hd_blue_clean    = np.copy(fluxes_hd_blue)

rcParams['figure.figsize'] = 12, 7 # Make plot small

for i in range(0,len(fluxes_pluto_red)):
  fluxes_pluto_red_clean[i] = clean_artifacts(fluxes_pluto_red[i])

for i in range(0,len(fluxes_pluto_blue)):
  fluxes_pluto_blue_clean[i] = clean_artifacts(fluxes_pluto_blue[i])

for i in range(0,len(fluxes_hd_red)):
  fluxes_hd_red_clean[i] = clean_artifacts(fluxes_hd_red[i])

for i in range(0,len(fluxes_hd_blue)):
  fluxes_hd_blue_clean[i] = clean_artifacts(fluxes_hd_blue[i])


##########
# Now merge and resample all of the spectra to a common baseine
##########

evenings_spect_pluto = sort(list(set(date_evening[w_spect_pluto_blue])))
evenings_spect_hd    = sort(list(set(date_evening[w_spect_hd_blue])))

num_spectra_pluto = size(evenings_spect_pluto) 	# Get list of unique nights on which we took Pluto spectra
num_spectra_hd    = size(evenings_spect_hd) 	# Get list of unique nights on which we took HD spectra

#   num_spectra_pluto = len(fluxes_pluto_blue) / 2
#   num_spectra_hd    = len(fluxes_hd_blue) / 2

fluxes_pluto_merged = []
fluxes_hd_merged = []

# Create the output wavelength array. 

DO_FORCE_FULL_RANGE = True

if (DO_FORCE_FULL_RANGE):
  wavelength_merged = np.arange(3000, 9500, dw)
else:
  wavelength_merged = np.arange(np.trunc(amin(wavelengths_pluto_blue)), np.trunc(amax(wavelengths_pluto_red)), dw)

print "Now creating output spectra (" + repr(np.amin(wavelength_merged)) + " .. " + repr(np.amax(wavelength_merged)) + \
  " AA) at %.2f A; originals are at %.2f A" % (dw, wavelengths_hd_red[0][1] - wavelengths_hd_red[0][0])

# NB: This %-style formatting is deprecated - should use .format instead.
 
#   for i in range(num_spectra_pluto): # Loop over the nights
#       print "Merging Pluto spectrum " + repr(i) + "/" + repr(num_spectra_pluto)
#       f = merge_four_rss(
#                   fluxes_pluto_blue_clean[2*i],   wavelengths_pluto_blue[2*i],
#                   fluxes_pluto_red_clean[2*i],    wavelengths_pluto_red[2*i],
#                   fluxes_pluto_blue_clean[2*i+1], wavelengths_pluto_blue[2*i+1],
#                   fluxes_pluto_red_clean[2*i+1],  wavelengths_pluto_red[2*i+1],
#                   wavelength_merged)               
#       fluxes_pluto_merged.append(f)

# Create variables to test multiple_rss

#       fluxes_blue_in      =  np.array([fluxes_pluto_blue_clean[2*i], fluxes_pluto_blue_clean[2*i+1], fluxes_pluto_blue_clean[2*i+3]])
#       fluxes_red_in       =  np.array([fluxes_pluto_red_clean[2*i], fluxes_pluto_red_clean[2*i+1]])
#       wavelengths_blue_in =  np.array([wavelengths_pluto_blue[2*i], wavelengths_pluto_blue[2*i+1], wavelengths_pluto_blue[2*i+1]])
#       wavelengths_red_in  =  np.array([wavelengths_pluto_red[2*i], wavelengths_pluto_red[2*i+1]])
#   
#       num_red_in  = shape(fluxes_red_in)[0]		# Number of red spectra passed
#       num_blue_in = shape(fluxes_blue_in)[0]		# Number of blue spectra passed
#       length_spectrum = shape(fluxes_blue_in)[1]
#       length_spectrum_merged = shape(wavelength_merged)[0]
#   #   
#   #   # Create the output where the drizzled fluxes go
#   #   
#       fluxes_red      = np.zeros((num_red_in,  length_spectrum_merged), dtype = float)
#       fluxes_blue     = np.zeros((num_blue_in, length_spectrum_merged), dtype = float)
#       wavelengths_blue = np.zeros((num_blue_in, length_spectrum), dtype = float)
#       wavelengths_red  = np.zeros((num_red_in,  length_spectrum), dtype = float)
#   
#       for i in range(num_red_in):
#         print 'Drizzling red ' + repr(i)
#         fluxes_red[i,:] = drizzle_flux_wavelength(fluxes_red_in[i,:], wavelengths_red_in[i,:], wavelength_merged)
#   
#       for i in range(num_blue_in):
#         print 'Drizzling blue ' + repr(i)
#         fluxes_blue[i,:] = drizzle_flux_wavelength(fluxes_blue_in[i,:], wavelengths_blue_in[i,:], wavelength_merged)

# Now assemble an array with all of the reds (or blues) for that night
# fluxes_pluto_red has the data to be assembled (30 arrays)
# w_spect_pluto_red is a 30-element array that lists -- in order -- the indices in the large array, for the files that are in the fluxes_pluto_red array

# Get a proper list of the spectra that are from a given evening, and are {Pluto | HD}, and {Red | Blue}, and otherwise in our 'good' list
# Problem: this 'fluxes_pluto_blue_clean

for e in evenings_spect_pluto:

# Set up a list for each spectrum

      fluxes_pluto_blue_e = []
      fluxes_pluto_red_e  = []
      wavelengths_pluto_blue_e = []
      wavelengths_pluto_red_e  = []

#    	Loop over all of the spectra in total. For each one, see if its index is on our 'good' list, and it's for tonight.
#      	If so, add that spectrum itself to a list. 
#       There certainly is a way to do this with array indexing, but after an hour of staring at it, I gave up and just used a 
#       brute-force loop.

      for i in range(w_spect_pluto_blue.size):
        if (date_evening[w_spect_pluto_blue[i]] == e):
	  fluxes_pluto_blue_e.append(fluxes_pluto_blue_clean[i])
	  wavelengths_pluto_blue_e.append(wavelengths_pluto_blue[i])

      for i in range(w_spect_pluto_red.size):
        if (date_evening[w_spect_pluto_red[i]] == e):
	  fluxes_pluto_red_e.append(fluxes_pluto_red_clean[i])
	  wavelengths_pluto_red_e.append(wavelengths_pluto_red[i])

      print "For evening of " + e + ": Pluto = "    + repr(shape(fluxes_pluto_red_e)[0])    + " red + " + repr(shape(fluxes_pluto_blue_e)[0]) + " blue; " 

      wavelengths_pluto_red_e = np.array(wavelengths_pluto_red_e)
      wavelengths_pluto_blue_e = np.array(wavelengths_pluto_blue_e)

      fluxes_pluto_blue_e = np.array(fluxes_pluto_blue_e)
      fluxes_pluto_red_e = np.array(fluxes_pluto_red_e)

#         w_fluxes_red_e_pluto  = set(where(date_evening == e)[0]).intersection(set(w_spect_pluto_red[:]))
#         w_fluxes_blue_e_pluto = set(where(date_evening == e)[0]).intersection(set(w_spect_pluto_blue[:]))
#   
#         fluxes_pluto_blue_e =      np.array(fluxes_pluto_blue_clean[list(w_fluxes_blue_e_pluto)])
#         wavelengths_pluto_blue_e = np.array(np.array(wavelengths_pluto_blue) [list(w_fluxes_blue_e_pluto)])
#         fluxes_pluto_red_e =       np.array(fluxes_pluto_red_clean [list(w_fluxes_red_e_pluto)])
#         wavelengths_pluto_red_e =  np.array(np.array(wavelengths_pluto_red)  [list(w_fluxes_red_e_pluto)])

      f = merge_multiple_rss(			# Merge the Pluto spectra from a given night
                fluxes_pluto_blue_e,
		wavelengths_pluto_blue_e, 
                fluxes_pluto_red_e,
		wavelengths_pluto_red_e,
		wavelength_merged)

      fluxes_pluto_merged.append(f)

      print "Just finished creating merged Pluto spectrum for " + e
      print


for e in evenings_spect_hd:

      fluxes_hd_blue_e    = []
      fluxes_hd_red_e     = []
      wavelengths_hd_blue_e    = []
      wavelengths_hd_red_e     = []

      for i in range(w_spect_hd_blue.size):
        if (date_evening[w_spect_hd_blue[i]] == e):
	  fluxes_hd_blue_e.append(fluxes_hd_blue_clean[i])
	  wavelengths_hd_blue_e.append(wavelengths_hd_blue[i])

      for i in range(w_spect_hd_red.size):
        if (date_evening[w_spect_hd_red[i]] == e):
	  fluxes_hd_red_e.append(fluxes_hd_red_clean[i])
	  wavelengths_hd_red_e.append(wavelengths_hd_red[i])

      print "For evening of " + e + ": HD = " + repr(shape(fluxes_hd_red_e)[0]) + " red + " + repr(shape(fluxes_hd_blue_e)[0]) + " blue."

      wavelengths_hd_red_e = np.array(wavelengths_hd_red_e)
      wavelengths_hd_blue_e = np.array(wavelengths_hd_blue_e)

      fluxes_hd_blue_e = np.array(fluxes_hd_blue_e)
      fluxes_hd_red_e = np.array(fluxes_hd_red_e)

      f = merge_multiple_rss(			# Merge the HD spectra from a given night
                fluxes_hd_blue_e, 
                wavelengths_hd_blue_e, 
                fluxes_hd_red_e, 
                wavelengths_hd_red_e, 
		wavelength_merged)

      fluxes_hd_merged.append(f)

      print "Just finished HD spectrum for " + e
      print

#         quit

#   for i in range(num_spectra_hd):  # Should be 5 of these
#       print "Merging HD spectrum " + repr(i) + "/" + repr(num_spectra_hd)
#       f = merge_four_rss(
#                   fluxes_hd_blue_clean[2*i],   wavelengths_hd_blue[2*i],
#                   fluxes_hd_red_clean[2*i],    wavelengths_hd_red[2*i],
#                   fluxes_hd_blue_clean[2*i+1], wavelengths_hd_blue[2*i+1],
#                   fluxes_hd_red_clean[2*i+1],  wavelengths_hd_red[2*i+1],
#                   wavelength_merged)               

# stop

# Clean the regions around the band gaps. Do this by removing gaps themselves, 
# as well as # any data within a few pixels of them. 
# ** No, I'd really rather not do this... it will make plots prettier, but will lose the CH4 data that I need.

# is_zero = (fluxes_pluto_merged[0] == 0)
# is_zero_wider = (is_zero) | (np.roll(is_zero,5)) | (np.roll(is_zero,-5))
# not_zero_wider = is_zero_wider == False
# 
# bins_good = not_zero_wider  # Boolean array, one flag per wavelength bin, as to whether to plot or not


# Print to the screen the location of the chip gaps. Note that my reduction makes these a bit wider than they really are.

# bins_good_edge_upper = (bins_good == True) & (np.roll(bins_good, 1) == False)
# bins_good_edge_lower = (bins_good == True) & (np.roll(bins_good,-1) == False)
# 
# wavelength_edge_upper = (wavelength_merged[bins_good_edge_upper])[1:]
# wavelength_edge_lower = (wavelength_merged[bins_good_edge_lower])[:-1]
# 
# print "Gap #  Lower  Upper"
# for i in range(wavelength_edge_lower.size):
#     print i, wavelength_edge_lower[i], wavelength_edge_upper[i]

# Now make a plot of all the Pluto spectra

rcParams['figure.figsize'] = 20, 10 # Make plot normal size

if (year_obs == 2014):
  index_good_spectra_pluto = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
  ylim = (0, 2)

if (year_obs == 2015):
  index_good_spectra_pluto = {0, 1, 2, 3, 4, 5, 6, 7} 
  ylim = (0, 1.5)
  
for i in range(num_spectra_pluto):
  bins_good = (fluxes_pluto_merged[i] > 0)					# Skip a few spectral bins
  if i in index_good_spectra_pluto:
      plot(wavelength_merged[bins_good], fluxes_pluto_merged[i][bins_good] + i*0.1, color=colors[i])
      plt.text(amax(wavelength_merged), 
               np.median(fluxes_pluto_merged[i][bins_good][-100:-1]) + i*0.1, 
#                  '  ' + date_evening[w_spect_pluto_blue[i*2]]
             '  ' + evenings_spect_pluto[i])					# Or is this the right date?

plt.title('Pluto, SALT ' + repr(year_obs), fontsize=24)
plt.xlim((3400, 9700))
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim(ylim)
file_out = dir_data + '/spectra_pluto.png'
plt.savefig(file_out)
print "Wrote: " + file_out
show()

##########
# Make a plot of all the HD spectra
##########

rcParams['figure.figsize'] = 20, 10 # Make plot normal size

if (year_obs == 2014):
  index_good_spectra_hd = {0,1,2,3, 4}
  ylim = (0, 1.1)

if (year_obs == 2015):
  index_good_spectra_hd = {0,1,2,3,4,5,6,7,8}	# 2015: Write all HD spectra
  ylim = (0, 1.5)

for i in range(num_spectra_hd):
  if i in index_good_spectra_hd:  
    plot(wavelength_merged[bins_good], fluxes_hd_merged[i][bins_good] + i*0.1, color=colors[i])
    plt.text(amax(wavelength_merged), 
             np.median(fluxes_hd_merged[i][bins_good][-100:-1]) + i*0.1, 
#                '  ' + date_evening[w_spect_hd_blue[i*2]])			# Is this the right date?
             '  ' + evenings_spect_hd[i])					# Or is this the right date?
plt.title('HD 146233, SALT ' + repr(year_obs), fontsize=24)
plt.ylim(ylim)
plt.xlim((3400, 9700))
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
file_out = dir_data + '/spectra_hd.png'
plt.savefig(file_out)
print "Wrote: " + file_out
show()

# Make a plot of the ratios

rcParams['figure.figsize'] = 20, 10 

if (year_obs == 2014):
  index_good_spectra_hd = {0,1,2,3, 4}
  ylim = (0.35,4)
#   index_good_spectra_pluto = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11} # 5 is messed up; 11 is noisy

if (year_obs == 2015):
  index_good_spectra_hd = {0,1,2,3, 4, 5, 6, 7}
  index_good_spectra_pluto = {0, 1, 2, 3, 4, 5, 6, 7, 8}
  ylim = (0.35,3)

for i in range(num_spectra_pluto):
  if i in index_good_spectra_pluto:
      plot(wavelength_merged[bins_good], fluxes_pluto_merged[i][bins_good] / fluxes_hd_merged[2][bins_good] + 
        i*0.25, color=colors[i])
      plt.text(amax(wavelength_merged), 
             np.median((fluxes_pluto_merged[i] / fluxes_hd_merged[2])[bins_good][-100:-1]) + i*0.25, 
#                '  ' + date_evening[w_spect_pluto_blue[i*2]])
             '  ' + evenings_spect_pluto[i])
plt.title('Pluto / HD 146233, SALT ' + repr(year_obs), fontsize=24)
plt.ylim(ylim)
plt.xlim((3400, 9700))
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)

file_out = dir_data + '/spectra_ratios.png'
plt.savefig(file_out)
print "Wrote: " + file_out
show()

#   quit

##########
# Write text file for output
# Be sure to write the filename using the date of the *evening* of the obs, not the proper UT date itself.
##########

DO_WRITE_TXT = True

if (DO_WRITE_TXT):

    for i in range(num_spectra_hd):
	d = np.zeros((2,wavelength_merged.size))
	d[0,:] = wavelength_merged
	d[1,:] = fluxes_hd_merged[i]
	d = d.transpose()
	date_str = evenings_spect_hd[i]
	file = dir_data + '/' + 'spect_hd_merged_' + date_str + '.txt'
	savetxt(file, d, delimiter=', ', header = 'Wavelength_[angstroms], Intensity_[arbitrary]')
	print 'Wrote: ' + repr(i) + ': ' + file
	
    for i in range(num_spectra_pluto):
	d = np.zeros((2,wavelength_merged.size))
	d[0,:] = wavelength_merged
	d[1,:] = fluxes_pluto_merged[i]
	d = d.transpose()
	date_str = evenings_spect_pluto[i]
	file = dir_data + '/' + 'spect_pluto_merged_' + date_str + '.txt'
	savetxt(file, d, delimiter=', ', header = 'Wavelength_[angstroms], Intensity_[arbitrary]')
	print 'Wrote: ' + repr(i) + ': ' + file
	
#####################

quit

DO_OTHER = False

if (DO_OTHER):
    i = 8

    plot(wavelengths_pluto_blue[i], fluxes_pluto_blue[i] / amax(fluxes_pluto_blue[i]) + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
    plot(wavelengths_pluto_red[i], fluxes_pluto_red[i]   / amax(fluxes_pluto_red[i])  + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
    plt.text(9400, i-i%2 + 0.3, date_evening[w_spect_pluto_blue[i]])

    i = 9

    plot(wavelengths_pluto_blue[i], fluxes_pluto_blue[i] / amax(fluxes_pluto_blue[i]) + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
    plot(wavelengths_pluto_red[i], fluxes_pluto_red[i]   / amax(fluxes_pluto_red[i])  + (i-i%2) + 0.1*(i%2), color=colors[i % len(colors)])
    plt.text(9400, i-i%2 + 0.3, date_evening[w_spect_pluto_blue[i]])

    plt.xlim((7000,8000))
#       show()


