#! /usr/bin/env python

"""
.. class:: HandleReducedELGPlate

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class handleReducedPlate is dedicated to handling reduced ELG plates from eBOSS : measuring line fluxes and assigning redshift flags.

"""
import numpy as n
import astropy.io.fits as fits
from scipy.interpolate import interp1d
import os
from os.path import join

run2d = "v5_9_1" # os.environ['RUN2D']
run1d = "v5_9_1" # os.environ['RUN1D']
		
class HandleReducedELGPlate:
	"""
	Now outdated by the SDSS4 svn product : elgredshiftflag
	
	Loads the environement proper to the SDSS survey :

        :param plate: plate number
        :param mjd: modified julian date
	
	"""
	def __init__(self,plate = 8123, mjd = 56931):
		self.plate = plate
		self.mjd = mjd
		
	def loadPlate(self):
		"""
		Opens the plate files: spPlate, spZbest. In the case one isworking on the Utah cluster.
		"""
		spfile = join( os.environ['BOSS_SPECTRO_REDUX'] , run2d , str(self.plate) , "spPlate-"+ str(self.plate) +"-"+ str(self.mjd) +".fits" )
		zbfile = join( os.environ['BOSS_SPECTRO_REDUX'] , run2d , str(self.plate) , run1d , "spZbest-" + str(self.plate) +"-"+ str(self.mjd) +".fits" )
		self.outputFile = join(os.environ['BOSS_SPECTRO_REDUX'] , run2d , str(self.plate) , run1d ,"spZ_ELGflag-" + str(self.plate) +"-"+ str(self.mjd) +".fits")
		# join(os.environ['BOSS_SPECTRO_REDUX'], run2d, str(self.plate), run1d, "spZ_ELGflag-" + str(self.plate) +"-"+ str(self.mjd) +".fits")
		# opens spPlate file
		hdulist = fits.open(spfile)
		c0 = hdulist[0].header['coeff0']
		c1 = hdulist[0].header['coeff1']
		npix = hdulist[0].header['naxis1']
		self.wavelength = 10.**(c0 + c1 * n.arange(npix))
		self.flux = hdulist[0].data
		self.fluxErr = hdulist[1].data**(-0.5)
		self.goodPix = (hdulist[1].data>0)
		print self.flux, self.goodPix
		hdulist.close()
		# opens spZbest file
		hdulist = fits.open(zbfile)
		self.zstruc = hdulist[1].data
		hdulist.close()
		hdulist = 0
		# defines what are galaxies 
		self.selection = (self.zstruc['Z']>0.) & (self.zstruc['Z'] > self.zstruc['Z_ERR'])
		self.Ngalaxies = len((self.selection).nonzero()[0])
		print "data loaded, Ngalaxy=", self.Ngalaxies

	def loadSpec(self,fiber):
		"""
		Opens the plate files: spPlate, spZbest. In the case one isworking on the Utah cluster.
		"""
		spfile = join( os.environ['BOSS_SPECTRO_REDUX'] , run2d , str(self.plate) , "spPlate-"+ str(self.plate) +"-"+ str(self.mjd) +".fits" )
		print spfile
		# opens spPlate file
		hdulist = fits.open(spfile)
		c0 = hdulist[0].header['coeff0']
		c1 = hdulist[0].header['coeff1']
		npix = hdulist[0].header['naxis1']
		goodPix = (hdulist[1].data[fiber-1]>0)
		self.wavelength = 10.**(c0 + c1 * n.arange(npix))[goodPix]
		self.flux = hdulist[0].data[fiber-1][goodPix]
		self.fluxErr = hdulist[1].data[fiber-1][goodPix]**(-0.5)
		hdulist.close()
		# opens spZbest file

	def save_result(self):
		"""
		Saves the results into a spZ_ELG-file.fits
		"""
		prihdr = fits.Header()
		prihdr['PLATE'] = self.plate
		prihdr['MJD'] = self.mjd
		prihdu = fits.PrimaryHDU(header=prihdr)

		thdulist = fits.HDUList([prihdu, self.tbhdu])
		os.system('rm '+ self.outputFile )
		thdulist.writeto( self.outputFile )
