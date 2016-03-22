# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:27:56 2015

@author: throop
"""

# Idea: this is a general class for processing SALT data.
# The pipeline is basically done, so this class is used for data *analysis*, not
# processing.
#
# Sample usage:
#

# a = HBTPlutoSalt('2015-04-03')
# (flux,wavelength) = HBTPlutoSalt.readspect()

import numpy as np

class HBTSaltPluto:

# Sample filename: 'spect_pluto_merged_2015-05-03.txt'. Pass this entire thing to routine.
    
    def __init__(self, filename):
        self.telescope = 'SALT'
        self.country = 'South Africa'
        self.filename = filename
        self.path_template = '/Users/throop/Data/SALT_Pluto_YEAR/product/'
        if self.filename.find('2015') > 0: # .find() returns -1 if not found, which is boolean TRUE...
            self.year = 2015
        if filename.find('2014') > 0:
            self.year = 2014
        
        self.pathname = self.path_template.replace('YEAR', repr(self.year))
                
    def readspect(self):

# Read a processed SALT spectrum from a .txt file, and return 
# a tuple of NP arrays: (wavelength, flux) [in order of X,Y like plotting]  
      
        d = np.loadtxt(self.pathname + self.filename, delimiter=',')

        flux = []
        wavelength = []
        
        flux.append(d[:,1])
        wavelength.append(d[:,0])

        flux = np.array(flux)
        wavelength = np.array(wavelength)
        
        return (wavelength, flux)
            
    def area(self):
        return 3.14 * (self.diameter/2)**2
        
    def d2(self):
        return self.diameter
