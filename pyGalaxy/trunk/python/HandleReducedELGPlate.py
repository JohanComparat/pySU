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
from lineListAir import *
import LineFittingLibrary as lineFit
from os.path import join

class HandleReducedELGPlate:
	"""
	Now outdated by the SDSS4 svn product : elgredshiftflag
	
	Loads the environement proper to the SDSS survey :

        :param plate: plate number
        :param mjd: modified julian date
	
	"""
	def __init__(self,plate = 8123, mjd = 56931, dV=-9999.99, fitWidth = 40. ):
		self.plate = plate
		self.mjd = mjd
		self.dV = dV
		self.fitWidth = fitWidth
		self.lfit  =  lineFit.LineFittingLibrary(fitWidth = self.fitWidth)

	def loadPlate(self):
		"""
		Opens the plate files: spPlate, spZbest. In the case one isworking on the Utah cluster.
		"""
		spfile = join( os.environ['BOSS_SPECTRO_REDUX'] , os.environ['RUN2D'] , str(self.plate) , "spPlate-"+ str(self.plate) +"-"+ str(self.mjd) +".fits" )
		zbfile = join( os.environ['BOSS_SPECTRO_REDUX'] , os.environ['RUN2D'] , str(self.plate) , os.environ['RUN1D'] , "spZbest-" + str(self.plate) +"-"+ str(self.mjd) +".fits" )
		self.outputFile = join(os.environ['BOSS_SPECTRO_REDUX'] , os.environ['RUN2D'] , str(self.plate) , os.environ['RUN1D'] ,"spZ_ELGflag-" + str(self.plate) +"-"+ str(self.mjd) +".fits")
		# join(os.environ['BOSS_SPECTRO_REDUX'], os.environ['RUN2D'], str(self.plate), os.environ['RUN1D'], "spZ_ELGflag-" + str(self.plate) +"-"+ str(self.mjd) +".fits")
		# opens spPlate file
		hdulist = fits.open(spfile)
		c0 = hdulist[0].header['coeff0']
		c1 = hdulist[0].header['coeff1']
		npix = hdulist[0].header['naxis1']
		self.wavelength = 10.**(c0 + c1 * n.arange(npix))
		self.flux = hdulist[0].data
		self.fluxErr = hdulist[1].data**(-0.5)
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

	def fitLineFluxes(self):
		"""
		Fits the line fluxes on each spectra.
		"""
		output = n.ones_like(n.empty([self.Ngalaxies,195]))*self.dV
		for jj,ii in enumerate(n.arange(1000)[self.selection]) :
			redshift = self.zstruc['Z'][ii]
			fiber = self.zstruc['FIBERID'][ii]
			flA,flErrA=self.flux[ii],self.fluxErr[ii]
			ok=(n.isnan(flA)==False)&(n.isnan(flErrA)==False)
			wl,fl,flErr=self.wavelength[ok],flA[ok]*1e-17,flErrA[ok]*1e-17
			d_out,m,h=[],[],[]
			# fits [OII] doublet
			datI,mI,hI=self.lfit.fit_Line_OIIdoublet(wl,fl,flErr, a0= n.array([ float(doubletList[0]), float(doubletList[2]) ]) * (1+ redshift ),lineName="O2_3728", p0_sigma=0.1, model = "gaussian")
			d_out.append(datI)
			m.append(mI)
			h.append(hI)
			# now fits the emission lines
			for li in lineList:
				datI,mI,hI=self.lfit.fit_Line(wl,fl,flErr, float(li[0])*(1+redshift), lineName=li[1], continuumSide=li[2], model="gaussian", p0_sigma=1)
				d_out.append(datI)
				m.append(mI)
				h.append(hI)
			# now fits the recombination lines 
			for li in recLineList:
				datI,mI,hI=self.lfit.fit_Line(wl,fl,flErr, float(li[0])*(1+redshift), lineName=li[1], continuumSide=li[2],model="gaussian")
				d_out.append(datI)
				m.append(mI)
				h.append(hI)
			# now outputs the result table
			heading="".join(h)
			out=n.hstack((d_out))
			out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*self.dV
			output[jj]=out

		# creates the fits table containing the fits information
		self.columns = [fits.Column(name="FIBERID",format="I", array = self.zstruc['FIBERID'][self.selection]) ]
		#self.new_columns = tbhdu.data.columns
		# adds all the line fits columns
		colNames = heading.split()
		for kk in range(len(colNames)):
			self.columns.append( fits.Column(name=colNames[kk],format='D', array= output.T[kk] )  )
			#self.new_columns += fits.Column(name=colNames[ii],format='D', array= output.T[ii] )

		self.tbhdu = fits.BinTableHDU.from_columns(self.columns)

	def assign_ELG_redshift_flag(self,SNRhigh=5., SNRlow=2.):
		"""
		Assigns redshift flags to ELGs by fitting emission lines.

		"""
		lines=n.array(["O2_3728", "Ne3_3869", "O3_4960", "O3_5007", "N2_6549", "N2_6585", "S2_6718", "S2_6732", "H1_1216", "H1_3970", "H1_4102", "H1_4341", "H1_4862", "H1_6564"])
		# constructs the high SNR and low SNR detection arrays
		lMin=3600.
		lMax=10350.
		detHArr,detLArr,lineWv=[],[],[]
		for line in lines:
			flux=self.tbhdu.data[line+'_flux']
			flux_err=self.tbhdu.data[line+'_fluxErr']
			if line=="O2_3728":
				lw=(self.tbhdu.data[line+"_a0a"]+self.tbhdu.data[line+"_a0b"])/2.
			else:
				lw=self.tbhdu.data[line + '_a0']
			lwObs=(lw>lMin)&(lw<lMax)
			lineWv.append(lw)
			detH=(lwObs)&(flux>0)&(flux_err>0)&(flux>=SNRhigh * flux_err)
			detL=(lwObs)&(flux>0)&(flux_err>0)&(flux>SNRlow * flux_err)& (flux<SNRhigh*flux_err)
			detHArr.append(detH)
			detLArr.append(detL)

		detHArr=n.transpose(detHArr)
		detLArr=n.transpose(detLArr)
		lineWv=n.transpose(lineWv)
		# now constructs the line quality array
		zQ,lname=[],[]
		for i in range(len(detHArr)):
			if self.zstruc['CLASS'][self.selection][i]=="STAR":
				zQ.append(-1.)
				lname.append(-2)

			elif self.zstruc['Z_ERR'][self.selection][i]>0.005* (1+ self.zstruc['Z'][self.selection][i] ) or self.zstruc['ZWARNING'][self.selection][i] > 4 :
				zQ.append(-2.)
				lname.append(-3)

			elif len(detHArr[i].nonzero()[0])>=3 :
				zQ.append(4.5)
				lname.append(lines[detHArr[i]])

			elif len(detHArr[i].nonzero()[0])==2 :
				zQ.append(4.) # two lines
				lname.append(lines[detHArr[i]])

			elif len(detHArr[i].nonzero()[0])==1 and detHArr[i][0]==True :
				zQ.append(3.5) # case of the OII doublet only
				lname.append(lines[detHArr[i]])

			elif len(detHArr[i].nonzero()[0])==1 and len(detLArr[i].nonzero()[0])>=1 :
				zQ.append(3.)
				lname.append(n.hstack((lines[detHArr[i]],lines[detLArr[i]])))

			elif len(detHArr[i].nonzero()[0])==0 and len(detLArr[i].nonzero()[0])>=3 :
				zQ.append(2.5)
				lname.append(lines[detLArr[i]])

			elif len(detHArr[i].nonzero()[0])==1 :
				zQ.append(2.)
				lname.append(lines[detHArr[i]])

			elif len(detHArr[i].nonzero()[0])==0 and len(detLArr[i].nonzero()[0])==2  :
				zQ.append(1.5)
				lname.append(lines[detLArr[i]])

			elif len(detHArr[i].nonzero()[0])==0 and len(detLArr[i].nonzero()[0])==1 :
				zQ.append(1.)
				lname.append(lines[detLArr[i]])

			else :
				zQ.append(0.)
				lname.append(-9)

		zQ=n.array(zQ)

		self.tbhdu.columns += fits.Column(name="zQ",format='D', array= zQ )
		self.tbhdu = fits.BinTableHDU.from_columns(self.tbhdu.columns)
		#self.new_columns += fits.Column(name="zCont",format='D', array= zCont )

	def save_result(self):
		"""
		Saves the results into a spZ_ELG-file.fits
		"""

		# now writes the output
		# creates the hdu
		prihdr = fits.Header()
		prihdr['PLATE'] = self.plate
		prihdr['MJD'] = self.mjd
		prihdu = fits.PrimaryHDU(header=prihdr)

		thdulist = fits.HDUList([prihdu, self.tbhdu])
		os.system('rm '+ self.outputFile )
		thdulist.writeto( self.outputFile )
