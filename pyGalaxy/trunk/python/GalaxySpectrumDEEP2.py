"""
.. class:: GalaxySpectrumDEEP2

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySpectrumDEEP2 is dedicated to handling DEEP2 spectra

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
import glob
from scipy.optimize import curve_fit
from GalaxySurveyDEEP2 import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
from LineFittingLibrary import *
lfl = LineFittingLibrary()
from filterList import *

class GalaxySpectrumDEEP2:
    """
    Loads the environement proper to the DEEP2 survey.

    Two modes of operation : flux calibration or line fitting
            
    :param catalog_entry: an entry of the deep2 catalog
    :param survey: survey python class
    :param calibration: if the class is loaded with intention of flux calibrating the DEEP2 data.
    :param lineFits: if the class is loaded with intention of fitting line fluxes on the DEEP2 spectra."""
    def __init__(self,catalog_entry, survey=GalaxySurveyDEEP2(redshift_catalog="zcat.deep2.dr4.v3.LFcatalogTC.Planck15.fits",calibration = False) ,calibration=False,lineFits=True ):
        self.catalog_entry=catalog_entry
        self.mask=str(self.catalog_entry['MASK'])
        self.slit=self.catalog_entry['SLIT']
        self.objno=str(self.catalog_entry['OBJNO'])

        self.database_dir = os.environ['DATA_DIR']
        self.deep2_dir = join(self.database_dir,"DEEP2")
        self.deep2_catalog_dir = join(self.deep2_dir,"catalogs")
        self.deep2_spectra_dir = join(self.deep2_dir,"spectra")
        self.survey = survey

        if calibration :
            self.path_to_spectrum = glob.glob(join(self.deep2_spectra_dir , self.mask ,'*', '*' + self.objno + '*.fits'))

        if lineFits :
            self.path_to_spectrum = glob.glob(join(self.deep2_spectra_dir , self.mask, '*', '*' + self.objno + '*_fc_tc.dat'))

        #print "path to spectrum", self.path_to_spectrum, self.mask, self.objno

        
    def openObservedSpectrum(self):
        """Loads the observed spectrum in counts."""
        hdS=fits.open(self.path_to_spectrum[0])
        self.airmass = hdS[1].header['AIRMASS']
        # blue spectrum
        self.dB=hdS[1].data
        # red pectrum
        self.dR=hdS[2].data
        self.chipNO=(hdS[1].header['CHIPNO']-1)%4
        #print hdS[1].header['CHIPNO']-1, self.chipNO
        hdS.close()
        lb=n.hstack((self.dB['LAMBDA'][0],self.dR['LAMBDA'][0]))
        self.lambdSwitch=n.max(self.dB['LAMBDA'][0])
        self.pixSampled=n.arange(2*4096)[(lb>6000)&(lb<10000)]
        self.lambd=lb[(lb>6000)&(lb<10000)]
        self.spec=n.hstack((self.dB['SPEC'][0],self.dR['SPEC'][0]))[(lb>6000)& (lb<10000)]
        self.ivar=n.hstack((self.dB['IVAR'][0],self.dR['IVAR'][0]))[(lb>6000)& (lb<10000)]
        self.specErr=self.ivar**(-0.5)

    def openObservedSpectrumFC(self):
        """Loads the observed spectrum in counts.
        """
        self.wavelength,self.fluxl,self.fluxlErr = n.loadtxt(self.path_to_spectrum[0] , unpack=True )

    def correctQE(self):
        """Corrects from the quantum efficiency of the chips where the spectrum landed."""
        if self.lambd.max()-self.lambd.min() > 3000 or n.mean(self.lambd)<7300 or n.mean(self.lambd)>8300 :
            print "cannot QE correct"

        xravg = 8900
        yravg = 150
        correctionavg = self.survey.paramsEndr[0] + self.survey.paramsEndr[1] * self.lambd
        self.xavg = (self.lambd - xravg)/yravg 
        ok1 =  (self.xavg > 0) & ( self.xavg < 1)
        self.cor2avg = correctionavg*self.xavg + 1*(1-self.xavg)
        ok2=(ok1)&(self.cor2avg>1)
        self.cor2avg[(ok2==False)] = n.ones_like(self.cor2avg[(ok2==False)])

        #npixel=len(self.lambd)
        self.left=(self.lambd<=self.lambdSwitch) # n.arange(4096)
        self.right=(self.lambd>self.lambdSwitch) # n.arange(4096,4096*2,1)

        #xx_b=self.lambd[self.left]
        #xx_r=self.lambd[self.right]

        #corr_b = params[num,0] + params[num,1]*self.lambd[self.left] + params[num,2]*self.lambd[self.left]**2
        #corr_r = params[num+4,0] + params[num+4,1]*self.lambd[self.right] + params[num+4,2]*self.lambd[self.right]**2
        corr_b = 1./( self.survey.params.T[self.chipNO][0] + self.survey.params.T[self.chipNO][1] * self.lambd[self.left] + self.survey.params.T[self.chipNO][2]*self.lambd[self.left]**2 )
        corr_r = 1./( self.survey.params.T[self.chipNO+4][0] + self.survey.params.T[self.chipNO+4][1]* self.lambd[self.right] + self.survey.params.T[self.chipNO+4][2] *self.lambd[self.right]**2 )
        # print corr_b, corr_r, self.cor2avg
        # print "spectrum",self.spec

        self.specTMP=n.zeros_like(self.spec)
        self.specErrTMP=n.zeros_like(self.specErr)
        self.ivarTMP=n.zeros_like(self.ivar)

        self.specTMP[self.left]=self.spec[self.left]*corr_b
        self.specTMP[self.right]=self.spec[self.right]*corr_r* self.cor2avg[self.right]

        self.specErrTMP[self.left]=self.specErr[self.left]*corr_b
        self.specErrTMP[self.right]=self.specErr[self.right]*corr_r* self.cor2avg[self.right]

        self.ivarTMP[self.left]=self.ivar[self.left]/(corr_b*corr_b)
        self.ivarTMP[self.right]=self.ivar[self.right]/(corr_r*corr_r* self.cor2avg[self.right]*self.cor2avg[self.right] )

        self.specTMP=self.specTMP/self.survey.throughput.y[self.pixSampled]
        self.specErrTMP=self.specErrTMP/self.survey.throughput.y[self.pixSampled]
        self.ivarTMP=self.ivarTMP*self.survey.throughput.y[self.pixSampled]**2

    def correct_telluric_abs(self):
        """ Future function to correct the observed spectra from tellurica absorption bands. Not yet implemented. """
        correction = self.survey.telluric_A_band_fct(self.lambd)**(self.airmass**0.55) * self.survey.telluric_B_band_fct(self.lambd)**(self.airmass**0.55)
        self.specTMP=self.specTMP/correction
        self.specErrTMP=self.specErrTMP/correction
        self.ivarTMP=self.ivarTMP*correction**2

    def fluxCal(self):
        """Performs the flux calibration of the spectrum by converting counts to flux units with an interpolation between the B and the I band borad band photometry."""
        countr = n.sum(self.specTMP*self.survey.Rresponse(self.lambd))/ n.sum(self.survey.Rresponse( self.lambd))
        counti = n.sum(self.specTMP*self.survey.Iresponse(self.lambd))/n.sum( self.survey.Iresponse( self.lambd))
        # (in erg/s/cm^2/Hz)
        fluxr = 10**((self.catalog_entry['MAGR'] + 48.6)/(-2.5)) 
        fluxi = 10**((self.catalog_entry['MAGI'] + 48.6)/(-2.5))
        fpcr = fluxr / countr
        fpci = fluxi / counti
        effr = 6599.0889
        effi = 8135.4026
        x = [effr, effi]
        y = [fpcr, fpci]
        # print x, y
        if y[0]>0 and y[1]>0:
            pfits = curve_fit(self.survey.fun,n.log(x),n.log(y),p0=(-0.01,-68))
            fluxn_corr = n.e**( pfits[0][1] + n.log(self.lambd)*pfits[0][0] )
        elif y[0]>0 and y[1]<0:
            fluxn_corr=fpcr
        elif y[0]<0 and y[1]>0:
            fluxn_corr=fpci
        else :
            return "bad"

        self.fluxn = fluxn_corr * self.specTMP
        self.fluxnErr = fluxn_corr * self.specErrTMP
        self.ivar_fluxn=self.ivarTMP/fluxn_corr**2

        self.fluxl=self.fluxn*299792458.0 / (self.lambd**2 * 10**(-10))
        self.fluxlErr=self.fluxnErr *299792458.0 / (self.lambd**2 * 10**(-10))
        return fluxn_corr


    def writeFCspec(self):
        """Writes the flux-calibrated spectrum"""
        ff=open(self.path_to_spectrum[0][:-5]+"_fc_tc.dat",'w')
        n.savetxt(ff,n.transpose([self.lambd,self.fluxl,self.fluxlErr]))
        ff.close()

    def openCalibratedSpectrum(self):
        """Loads the flux calibrated spectrum in f lambda convention.
        """
        self.wavelength,self.fluxl,self.fluxlErr= n.loadtxt(self.path_to_spectrum[0], unpack=True)


	def plotFit(self, outputFigureNameRoot, ymin = 1e-19, ymax = 1e-17):
		"""
		Plots the spectrum and the line fits in a few figures
		"""
		ifl = lfl.flambda(self.catalog_entry['MAGI'], lambIcfht)
		ifl_max = lfl.flambda(self.catalog_entry['MAGI']+self.catalog_entry['MAGIERR'], lambIcfht)
		ifl_min = lfl.flambda(self.catalog_entry['MAGI']-self.catalog_entry['MAGIERR'], lambIcfht)

		rfl = lfl.flambda(self.catalog_entry['MAGR'], lambRcfht)
		rfl_max = lfl.flambda(self.catalog_entry['MAGR']+self.catalog_entry['MAGRERR'], lambRcfht)
		rfl_min = lfl.flambda(self.catalog_entry['MAGR']-self.catalog_entry['MAGRERR'], lambRcfht)
		
		p.figure(1,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr, linewidth=1, alpha= 0.4, label='spectrum')
		p.plot([lambIcfht,lambIcfht,lambIcfht],[ifl_min,ifl,ifl_max], 'r', label = 'magnitudes', lw=2)
		p.plot([lambRcfht, lambRcfht, lambRcfht], [rfl_min, rfl, rfl_max], 'r', lw=2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.savefig( outputFigureNameRoot + "-all.png" )
		p.clf()

		a0 = self.catalog_entry['O2_3728_a0']
		continu= self.catalog_entry['O2_3728_continu']
		aas =n.arange(self.catalog_entry['O2_3728_a0']-70, self.catalog_entry['O2_3728_a0']+70,0.1)
		flMod=lambda aa,sigma,F0,sh :continu+ lfl.gaussianLineNC(aa,sigma,(1-sh)*F0,a0)+lfl.gaussianLineNC(aa,sigma,sh*F0,a0)
		model = flMod(aas, self.catalog_entry['O2_3728_sigma'], self.catalog_entry['O2_3728_flux'], self.catalog_entry['O2_3728_share'] )# self.catalog_entry['O2_3728_share'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.axvline(self.catalog_entry['O2_3728_a0'],color='k', ls='dashed', label= 'obs')
		p.plot(aas, model,'g',label='model', lw=2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O2_3728_a0']-100, self.catalog_entry['O2_3728_a0']+100))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OII] 3727')
		p.savefig( outputFigureNameRoot + "-O2_3728.png")
		p.clf()

		a0 = self.catalog_entry['O3_5007_a0']
		continu= self.catalog_entry['O3_5007_continu']
		aas =n.arange(self.catalog_entry['O3_5007_a0']-70, self.catalog_entry['O3_5007_a0']+70,0.1)
		flMod=lambda aa,sigma,F0: lfl.gaussianLine(aa,sigma,F0,a0,continu)
		model = flMod(aas, self.catalog_entry['O3_5007_sigma'], self.catalog_entry['O3_5007_flux'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.axvline(self.catalog_entry['O3_5007_a0'],color='k', ls='dashed', label= 'obs')
		p.plot(aas, model,'g',label='model', lw =2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O3_5007_a0']-100, self.catalog_entry['O3_5007_a0']+100))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OIII] 5007')
		p.savefig( outputFigureNameRoot + "-O3_5007.png")
		p.clf()


