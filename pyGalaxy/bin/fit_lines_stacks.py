"""
This script fits emission lines on a stacked spectrum
"""
import os
import glob
import numpy as n
import astropy.io.fits as fits

from lineListAir import *
allLinesList = n.array([ [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"], [H1,H1_1216,"H1_1216","right"], [H1,H1_3970,"H1_3970","right"], [H1,H1_4102,"H1_4102","right"], [H1,H1_4341,"H1_4341","right"], [H1,H1_4862,"H1_4862","left"], [H1,H1_6564,"H1_6564","left"]])

doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])

# import the fitting routines
import LineFittingLibrary as lineFit


fileList=n.array(glob.glob("/home/comparat/database/*/products/emissionLineLuminosityFunctions/??_????/??_????-*-z*.fits"))

fileList2=n.array(glob.glob( "/home/comparat/database/*/products/emissionLineLuminosityFunctions/??_????/??_????-*-z*stack*+??.fits") )

stackSel = n.array([ el.find('stack') for el in fileList])

stackFiles = (stackSel != -1)
LFfiles = (stackSel == -1)

#fileList[stackFiles]
LFfileList = fileList[LFfiles]
stacksToUse = n.array([n.array(glob.glob(el[:-5]+"*stack*+??.fits")) for el in LFfileList ])

# Then a loop on the 2 lists 
# 
for ii in range(len(LFfileList)):
	if len(stacksToUse[ii])>0 :
		for jj in range(len(stacksToUse[ii])):
			mainFile = LFfileList[ii]
			stackFile = stacksToUse[ii][jj]
			print mainFile.split('/')[-1], stackFile.split('/')[-1]
			if stackFile.find('VVDS')>0 or stackFile.find('VIPERS')>0 :
				lfit  =  lineFit.LineFittingLibrary(fitWidth = 70.)
			if stackFile.find('DEEP2')>0 :
				lfit  =  lineFit.LineFittingLibrary(fitWidth = 40.)

			hdus = fits.open(stackFile)
			hdR = hdus[0].header
			hdu1 = hdus[1] # .data
			print "does the fits"
			wlA,flA,flErrA = hdu1.data['wavelength'], hdu1.data['meanWeightedStack'], hdu1.data['jackknifStackErrors']
			selection = (flA>0) & (hdu1.data['NspectraPerPixel']  > float( stackFile.split('_')[-5]) * 0.7 )
			wl,fl,flErr = wlA[selection], flA[selection], flErrA[selection] 

			data,h=[],[]

			fl = flA[selection]
			dat_mean,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl, flErr, a0= n.array([O2_3727,O2_3729]) , lineName="O2_3728", p0_sigma=1,model="gaussian")

			d_out=[]
			for kk in range(10):
				fl = hdu1.data['jackknifeSpectra'].T[kk][selection]
				d1,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl, flErr, a0= n.array([O2_3727,O2_3729]) , lineName="O2_3728", p0_sigma=1,model="gaussian")
				d_out.append(d1)

			d_out = n.array(d_out)
			err_out = n.std(d_out,axis=0)
			# assign error values :
			dat_mean[2] = err_out[2-1]
			dat_mean[4] = err_out[4-1]
			dat_mean[6] = err_out[6-1]
			data.append(dat_mean)
			h.append(hI)

			for li in allLinesList :
				# measure line properties from the mean weighted stack
				fl = flA[selection]
				dat_mean,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1], lineName=li[2], continuumSide=li[3], model="gaussian",p0_sigma=1)
				# measure its dispersion using the stacks
				d_out=[]
				for kk in range(len(hdu1.data['jackknifeSpectra'].T)):
					fl = hdu1.data['jackknifeSpectra'].T[kk][selection]
					d1,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1], lineName=li[2], continuumSide=li[3], model="gaussian",p0_sigma=1)
					d_out.append(d1)

				d_out = n.array(d_out)
				err_out = n.std(d_out,axis=0)
				# assign error values :
				dat_mean[2] = err_out[2-1]
				dat_mean[4] = err_out[4-1]
				dat_mean[6] = err_out[6-1]
				data.append(dat_mean)
				h.append(hI)

			heading="".join(h)
			out=n.hstack((data))
			out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*lfit.dV
			#output = n.array([ out ])
			#print "----------------", output.T[0], output.T[1], output
			colNames = heading.split()
			#col0 = fits.Column(name=colNames[0],format='D', array= output.T[0])
			#col1 = fits.Column(name=colNames[1],format='D', array= output.T[1])
			#cols = fits.ColDefs([col0, col1])
			#print colNames
			for ll in range(len(colNames)):
				hdR[colNames[ll]] = out.T[ll]
				#cols += fits.Column(name=colNames[ll],format='D', array= output.T[ll] )

			#tbhdu = fits.BinTableHDU.from_columns(cols)
			tbhdu = fits.BinTableHDU.from_columns( hdu1.columns )
			prihdu = fits.PrimaryHDU(header=hdR)
			thdulist = fits.HDUList([prihdu, tbhdu])
			os.system('rm '+ stackFile[:-5]+'_linesFitted.fits')
			thdulist.writeto(stackFile[:-5]+'_linesFitted.fits')


