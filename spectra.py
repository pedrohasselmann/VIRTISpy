# -*- coding: utf-8 -*-
#!/usr/bin/python

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

__version__="1.0.0"

__doc__='''
  ############################ MAIN PLOT TOOL ################################
  ### Author: Pedro Henrique A. Hasselmann
  ### Last Modified: March, 2016
  ###
  ### Imaging Spectrometer analysis tool
  ###
  ###
  ### Classes: Spectra
  ###
  ###
  ### Attributes:
  ###  
  ### 
  ###  
  ### dependencies:
  ###             numpy, scipy
  ###             
  ###
  #############################################################################
'''

# Global Import
from source import *

###########################
### SPECTROMETRIC IMAGE ###
###########################

class LowSpectra:
    '''
     Tool for calibrating and manipulating low-resolution Spectrometric images.
     
     Parameters
     ==========
     
     data         --> encoded VIRTIS DATA class (virtis.py)
     geo          --> encoded VIRTIS GEO class (virtis.py)
     solar_spec   --> 3-column array containing calibrated solar spectrum in the same spectral range. 
    '''
    import pandas as pd
    global pd


    def __init__(self, cubename, mission, cubegeoname=None, spectral_range=None,**args):   
       from source.spectra.virtis import read_VIRTIS
       from numpy import around
    
       self.data = read_VIRTIS(cubename, mission=mission)
    
       if cubegeoname != None: self.geo = read_VIRTIS(cubegeoname, mission=mission)
       
       self.cubename = cubename
       self.cubegeoname = cubegeoname
       
       if spectral_range:
         from numpy import where
         m1, m2 = spectral_range
         wv = self.data.getWave()
         i = where((wv >= m1) & (wv <= m2))[0]
         a1, a2 = i[0], i[-1]
         self.spindex = (a1, a2)
         print(a1, a2)
         self.data.cube = self.data.cube[a1:a2,...]
         self.data.table = self.data.table[a1:a2,...]
 
       self.wv_index = pd.Series(range(self.data.getWave().size),index=around(self.data.getWave()))


    def solar_calib(self, solar_spec):
       from math import pi
       R = float(self.data.header['SOLAR_DISTANCE'])/au_km #u.au.to(u.km)
       #R = 490196206.065/au_km
       print('Solar Distance: ',R)
       #self.solar = (pi*R**2)/solar_spec[:,-1]
       if hasattr(self, 'spindex'): solar_spec = solar_spec[self.spindex[0]:self.spindex[1]]
       self.data.cube = (pi*R**2)*(self.data.cube.T/solar_spec[:,-1]).T

       
    def polyfit(self, deg, window_div=2, weigth=None, norm_at=None, **args):
      '''
        1. We filter the spectra using Savitzky-Golay filter(scipy.savgol_filter)
        2. We trace a polynomial to compute the corresponding coefficients to each pixel.
        
        Visible spectral slope:
        S = 10^4 * (r_800 - r_450)/(r_450 *(800-450))
      '''
      import numpy.polynomial.polynomial as poly
      from numpy import float32, isnan, isinf, isneginf, inf, nan, median, where, diff, transpose
      from scipy import signal

      wv = self.data.getWave()
      x, y, z = self.data.cube.shape
      #print([x, y, z])
      data_table = self.data.cube.reshape(x,y*z)

      # Savitzky-Golay filter: Cleaning spectra
      refl_filtered = signal.savgol_filter(data_table, axis=0, window_length=int(wv.size/window_div)+1, polyorder=deg)


      # polynomial fit
      fit, res = poly.polyfit(wv, refl_filtered, deg=deg, full=True)
      line, res_line = poly.polyfit(wv, refl_filtered, deg=1, full=True)
      #self.refl_poly = poly.polyval(wv, fit).T.reshape(x,y,z)

      #import matplotlib.pyplot as plt
      #plt.plot(wv, data_table.mean(1), linewidth=3.5,label='Spectrum')
      #plt.plot(wv, refl_filtered.mean(1), linewidth=3.5,label='Filtered Savitzky-Golay')
      #plt.plot(wv, poly.polyval(wv, fit).T.mean(1), linewidth=3.5,label='poly-5')
      #plt.plot(wv, poly.polyval(wv, line).T.mean(1), linewidth=3.5,label='line')
      #plt.legend(loc=0)
      #plt.show()
     
      self.vis_slope = 1e4*(line[1,:]/poly.polyval(norm_at, line)).reshape(y,z)


      self.fit_poly = fit.reshape(deg+1,y,z)      
      self.refl_filtered = refl_filtered.reshape(x,y,z)
  
    # Photometric disk correction
    def diskcorr(diskfunc):
      assert self.geo
      pass
      
  
    # FITS Image by astropy
    def to_fits(self, parameter, label='spec_slope'):
      '''
         output data into FITS image.
      '''
      from astropy.io import fits
      from numpy import float32, isnan, isinf, isneginf, inf, nan
      
      image  = fits.PrimaryHDU(parameter.astype(float32))
      image.writeto(path.join(prod, label+'_'+self.cubename.split('/')[-1][:-4]+'.fit'), output_verify='warn')      


class HighSpectra:
    '''
     Tool for calibrating and manipulating high-resolution spectroscopic data.
     
     Parameters
     ==========
     
     data         --> encoded as 3-dim numpy.array
     solar_spec   --> 3-column array containing calibrated solar spectrum in the same spectral range. 
    '''

    def __init__(self, cubename, mission, cubegeoname=None, spectral_range=None,**args):   
       from source.spectra.virtis import read_VIRTIS
    
       self.data = read_VIRTIS(cubename, mission=mission)
    
       if cubegeoname != None: self.geo = read_VIRTIS(cubegeoname, mission=mission)
       
       self.cubename = cubename
       self.cubegeoname = cubegeoname
   
    def solar_calib(self, solar_spec):
       from math import pi
       R = float(self.data.header['SOLAR_DISTANCE'])/au_km #u.au.to(u.km)
       #R = 490196206.065/au_km
       print('Solar Distance: ',R)
       #self.solar = (pi*R**2)/solar_spec[:,-1]
       if hasattr(self, 'spindex'): solar_spec = solar_spec[self.spindex[0]:self.spindex[1]]
       self.data.cube = (pi*R**2)*(self.data.cube.T/solar_spec[:,-1]).T

       
    def polyfit(self, deg, window_div=2, weigth=None, **args):
      '''
        1. We filter the spectra using Savitzky-Golay filter(scipy.savgol_filter)
      '''
      import numpy.polynomial.polynomial as poly
      from numpy import float32, isnan, isinf, isneginf, inf, nan, median, where, diff
      from scipy import signal

      # Savitzky-Golay filter: Cleaning spectra
      refl_filtered = signal.savgol_filter(self.data.cube, axis=0, window_length=int(wv.size/window_div)+1, polyorder=deg)  
    
  
