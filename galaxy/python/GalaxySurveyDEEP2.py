"""
.. class:: GalaxySurveyDEEP2

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySurveyDEEP2 is dedicated to handling DEEP2 survey and the class GalaxySpectrumDEEP2 to handling its spectra.

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
from scipy.interpolate import interp1d
from MiscellanousFunctionsLibrary import *
import astropy.units as u

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

class GalaxySurveyDEEP2:
    """
    Loads the environement proper to the DEEP2 survey :
     * Defines all the proper paths in the database,
         * Opens the catalog and different calibration files,
         * Loads a list of the DEEP2 spectra.

    :param redshift_catalog: name of the DEEP2 redshift catalog (path to the fits file)
    :param calibration: if the class is loaded with intention of flux calibrating the DEEP2 data (boolean)"""
    def __init__(self,redshift_catalog="zcat.deep2.dr4.v4.fits", calibration=True, plots=True):
        self.redshift_catalog = redshift_catalog
        self.database_dir = '/data42s/comparat/firefly/v1_1_0'
        self.deep2_dir = join(self.database_dir,"DEEP2")
        self.deep2_catalog_dir = join(self.deep2_dir,"catalogs")
        self.deep2_spectra_dir = join(self.deep2_dir,"spectra")
        hd = fits.open(join(self.deep2_catalog_dir,self.redshift_catalog))
        self.catalog = hd[1].data
        hd.close()
        self.plots = plots

        if calibration==True :
            self.deep2_calib_dir = join(os.environ['GIT_PYSU'],"galaxy","data","Deep2_calib_files")
            self.paramsEndr = fits.open(join(self.deep2_calib_dir,"paramsendr.fits"))[0].data
            self.params = fits.open(join(self.deep2_calib_dir,"params.fits"))[0].data
            v0,v1 = n.loadtxt(join(self.deep2_calib_dir, "thr_go1200_80_og550.asc"), unpack = True) 
            self.throughput = interp1d( v0,v1 )

            self.telluric_A_band = fits.open(join(self.deep2_calib_dir,"aband.fits"))
            a_lambda = self.telluric_A_band[0].header['CRVAL1']+n.arange(self.telluric_A_band[0].header['NAXIS1'])*self.telluric_A_band[0].header['CDELT1']
            a_fluxn = self.telluric_A_band[0].data / ( ( n.sum(self.telluric_A_band[0].data[1400:1499])+n.sum(self.telluric_A_band[0].data[2500:2599]) )/200.)
            minab=1510
            maxab=2110
            errors=n.ones_like(a_lambda)
            errors[minab:maxab]=n.ones_like(a_lambda[minab:maxab])*1e20
            # fit a degree 3 polynomial to the band 
            a_coeff = n.polynomial.polynomial.polyfit( a_lambda, a_fluxn, deg = 1, w=1/errors)
            a_band_fit = n.polynomial.polynomial.polyval(a_lambda, a_coeff)
            a_flux = a_fluxn / a_band_fit
            self.telluric_A_band_fct = interp1d( n.hstack(([3000],a_lambda[minab:maxab],[10000])), n.hstack(([1.],a_flux[minab:maxab],[1.])) )
            if plots :
                p.figure(1,(5,5))
                p.plot(a_lambda,self.telluric_A_band[0].data,label='raw a Band data')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","a_band.pdf"))
                p.clf()

                p.figure(1,(5,5))
                p.plot(a_lambda,a_fluxn,label='raw a Band data normed')
                p.plot(a_lambda,a_band_fit,label='polynomial fit')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","a_band_normed.pdf"))
                p.clf()

                p.figure(1,(5,5))
                p.plot(a_lambda,a_flux,label='final a Band')
                p.plot(self.telluric_A_band_fct.x,self.telluric_A_band_fct.y,label='interpolation')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","a_band_final.pdf"))
                p.clf()

            self.telluric_B_band = fits.open(join(self.deep2_calib_dir,"bband.fits"))
            b_lambda = self.telluric_B_band[0].header['CRVAL1']+n.arange(self.telluric_B_band[0].header['NAXIS1'])*self.telluric_B_band[0].header['CDELT1']
            band1=(b_lambda<=6840)&(b_lambda>6840-25)
            band2=(b_lambda<=6960+25)&(b_lambda>6960)
            b_fluxn = self.telluric_B_band[0].data / ( ( n.sum(self.telluric_B_band[0].data[band1])+n.sum(self.telluric_B_band[0].data[band2]) )/200.)
            minab= (b_lambda>6650) # avoid the false feature (absorption spike) at 6556 A.
            bband_coverage = (b_lambda<=6960) & (b_lambda>6840) # the B band coverage
            lambda_for_fit = (minab) & (bband_coverage == False)
            errors=n.ones_like(b_lambda)
            errors[(lambda_for_fit == False)]=n.ones_like(b_lambda[(lambda_for_fit == False)])*1e20
            # fit a degree 3 polynomial to the band 
            b_coeff = n.polynomial.polynomial.polyfit( b_lambda, b_fluxn, deg = 1, w=1/errors)
            b_band_fit = n.polynomial.polynomial.polyval(b_lambda, b_coeff)
            b_flux = b_fluxn / b_band_fit
            self.telluric_B_band_fct = interp1d( n.hstack(([3000],b_lambda[bband_coverage],[10000])), n.hstack(([1.],b_flux[bband_coverage],[1.])) )
            if plots :

                p.figure(1,(5,5))
                p.plot(b_lambda,self.telluric_B_band[0].data,label='raw b Band data')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","b_band.pdf"))
                p.clf()

                p.figure(1,(5,5))
                p.plot(b_lambda,b_fluxn,label='raw b Band data normed')
                p.plot(b_lambda,b_band_fit,label='polynomial fit')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","b_band_normed.pdf"))
                p.clf()

                p.figure(1,(5,5))
                p.plot(b_lambda,b_flux,label='final b Band')
                p.plot(self.telluric_B_band_fct.x,self.telluric_B_band_fct.y,label='interpolation')
                p.xlabel('wavelength')
                p.ylabel('telluric band')
                p.legend(loc=3)
                p.savefig(join(self.deep2_spectra_dir,"plots","b_band_final.pdf"))
                p.clf()


            v0,v1 = n.loadtxt(join(self.deep2_calib_dir ,"Bresponse.txt"), unpack = True, usecols = (0,6))
            self.Bresponse = interp1d(v0,v1 )
            v0,v1 = n.loadtxt(join(self.deep2_calib_dir,"Rresponse.txt"), unpack = True, usecols = (0,6))
            self.Rresponse = interp1d(v0,v1 )	
            v0,v1 = n.loadtxt(join(self.deep2_calib_dir,"Iresponse.txt"), unpack = True, usecols = (0,6))
            self.Iresponse = interp1d(v0,v1 )
            self.fun = lambda x,a,b : a*x+b

        
    def computeLineLuminosity(self,line,distanceCorrection):
        """ computes the line luminosities for the line list given.
        :param catalog: fits catalog containing redshift, EBV and line fluxes
        :param line:
        """
        ebvCorrection=n.array([ 10**(0.4 *self.catalog['SFD_EBV'][i] * CalzettiLaw((1 + self.catalog['ZBEST'][i]) * line[1])) for i in range(len(self.catalog['ZBEST']))])
        flux=ebvCorrection*self.catalog[line[2]+'_flux']*u.erg/u.cm**2/u.s
        Luminosity=fits.Column(name=line[2]+"_luminosity",format="D", unit="erg/s", array=distanceCorrection*flux )
        LuminosityErr=fits.Column(name=line[2]+"_luminosityErr",format="D", unit="erg/s", array= self.catalog[line[2]+'_fluxErr']/ self.catalog[line[2]+'_flux']* distanceCorrection *flux)
        return Luminosity, LuminosityErr


