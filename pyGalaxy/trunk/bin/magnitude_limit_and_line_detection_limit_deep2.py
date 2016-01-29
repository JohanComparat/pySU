#! /usr/bin/env python

"""
This script fits emission lines on a stacked spectrum
"""
import os
import glob
import numpy as n
import astropy.io.fits as fits
import matplotlib.pyplot as p
# import the fitting routines
import LineFittingLibrary as lineFit
from scipy.interpolate import interp1d
import time
lfit  =  lineFit.LineFittingLibrary(fitWidth = 40.)

# CASE VVDS DEEP
iMags = n.arange(22,24.1,0.1)
snr5Flux = n.empty(len(iMags))
snr7Flux = n.empty(len(iMags))
snr5EW = n.empty(len(iMags))
snr7EW = n.empty(len(iMags))
for jj,mag in enumerate(iMags):
	print mag, time.time()
	signal_level =  lfit.flambda(mag,8300)
	# fixed noise level for the DEEP
	noise_level = lfit.flambda(24,8300) / 0.897 # 3.56*10**-19
	wavelength = n.arange(8000,8600,7)
	error = n.ones_like(wavelength) * noise_level
	flux_no_line = lfit.gaussianLine(wavelength,10,0,8300,signal_level)
	total = n.sum(flux_no_line)
	shares = n.arange(0.,1,0.05)
	data=n.empty((len(shares),12))
	for ii,share in enumerate(shares):
		print share, time.time()
		flux_with_line_clean = lfit.gaussianLine(wavelength,2,total*share,8300,total*(1-share)/len(wavelength))
		flux_with_line = flux_with_line_clean + noise_level * ( n.random.random(len(wavelength))*2.-1. )
		dat_mean,mI,hI=lfit.fit_Line(wavelength,flux_with_line,error,8300, lineName="O3_5007", continuumSide="left", model="gaussian",p0_sigma=1,DLC=100)
		data[ii]=dat_mean
		
	SNR =data.T[1]/data.T[2]
	snrBins = n.arange(0,100,2)
	snrVal = (snrBins[1:]+snrBins[:-1])/2.
	medianFlux = n.empty(len(snrVal))
	medianEW = n.empty(len(snrVal))
	for ii in range(len(snrVal)):
		sel = (snrBins[ii]<SNR)&(SNR<snrBins[ii+1])
		medianFlux[ii] = n.median(data.T[1][sel])
		medianEW[ii] = n.median(data.T[7][sel])

	snr_to_flux = interp1d(snrVal,medianFlux)
	snr_to_ew = interp1d(snrVal,medianEW)
	snr5Flux[jj] = snr_to_flux(5)
	snr7Flux[jj] = snr_to_flux(7)
	snr5EW[jj] = snr_to_ew(5)
	snr7EW[jj] = snr_to_ew(7)

n.savetxt("mag-limit-flux-ew-vvdsdeep.dat",n.transpose([iMags, snr5Flux, snr7Flux, snr5EW, snr7EW]))
n.savetxt("snr-flux-ew-relation-i24-vvdsdeep.dat",n.transpose([	snrVal, medianFlux, medianEW]))
	
p.plot(iMags,snr5Flux,'b+',label="SNR=5")
p.plot(iMags,snr7Flux,'r+',label="SNR=7")
p.xlabel('i mag limit')
p.ylabel('flux')
p.xlim((21.5,25))
p.ylim((1e-19,5e-16))
p.yscale('log')
p.grid()
p.legend(loc=3)
p.savefig("../../../../../snr-mag-flux-vvdsdeep.png")
p.clf()


p.plot(iMags,snr5EW,'b+',label="SNR=5")
p.plot(iMags,snr7EW,'r+',label="SNR=7")
p.xlabel('i mag limit')
p.ylabel('EW')
p.xlim((21.5,25))
#p.ylim((1e-19,5e-16))
p.yscale('log')
p.grid()
p.legend(loc=3)
p.savefig("../../../../../snr-mag-ew-vvdsdeep.png")
p.clf()

import sys
sys.exit()


from scipy.optimize import curve_fit
pl1 =lambda x, a0, a1 : a0 + a1*x
sel_fit=(SNR>5)&(SNR<100)
popt,pcov=curve_fit(pl1,n.log10(SNR[sel_fit]),n.log10(data.T[1][sel_fit]),p0=(-19,0.25))


p.plot(SNR,data.T[1],'b,')
p.plot(snrVal,medianFlux,'r--',lw=2,label='median value')
p.plot(snrVal[snrVal>5],10**pl1(n.log10(snrVal[snrVal>5]),popt[0],popt[1]),'k--',label=r'fit $10^{' +str(n.round(popt[0],2))+r'+'+ str(n.round(popt[1],2))+ r'(log_{10}(SNR))}$')
p.xlabel('SNR')
p.ylabel('flux')
p.xlim((1,100))
p.ylim((1e-19,5e-16))
p.yscale('log')
p.xscale('log')
p.grid()
p.legend(loc=2)
p.savefig("../../../../../snr-flux.png")
p.clf()
