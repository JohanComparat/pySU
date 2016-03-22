import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
from HandleSdssPlate import *
import glob	

stackDir = "/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG"
ggrid  = [21.8,22.5,22.8]
rzgrid = [0.0,0.8,1.0,2.0]
grgrid = [0.0,0.4,0.6,1.0]

def compareSpectrumToStack(entry, nameRoot="elg270_eboss17_", ggrid  = [21.8,22.5,22.8], rzgrid = [0.0,0.8,1.0,2.0], grgrid = [0.0,0.4,0.6,1.0]):
	# gets the spectrum
	ObsPlate = HandleReducedELGPlate(entry['PLATE'],entry['MJD'])
	ObsPlate.loadSpec(entry['FIBER'])
	# gets the stack
	if entry['index_g']>=0:
		suffix = "_g_"+str(n.round(ggrid[entry['index_g']],1))+"_rz_"+str(n.round(rzgrid[entry['index_rz']],1))+"_gr_"+str(n.round(grgrid[entry['index_gr']],1))
		stackName = join(stackDir, nameRoot + suffix + "_stack.fits")
		hdu = fits.open(stackName)[1].data
		sel = (hdu['NspectraPerPixel']>0.9*n.max(hdu['NspectraPerPixel']))
		# compares stakc and spectrum at REDSHIFT
		def getchi2(REDSHIFT):
			wlmin=n.min(hdu['wavelength'][sel]*(1+REDSHIFT))
			wlmax=n.max(hdu['wavelength'][sel]*(1+REDSHIFT))
			meanStack =interp1d(hdu['wavelength'][sel]*(1+REDSHIFT),hdu['medianStack'][sel])
			overlap=(ObsPlate.wavelength>wlmin)&(ObsPlate.wavelength<wlmax)
			x = ObsPlate.wavelength[overlap]
			y = ObsPlate.flux[overlap]
			yerr = ObsPlate.fluxErr[overlap]
			chi2mean = n.sum(((y - meanStack(x))/yerr)**2.)/len(x)
			return chi2mean
			
		return getchi2(entry['Z_1']), getchi2(entry['Z_2']), getchi2(entry['Z_3'])
	else:
		return -1,-1,-1


def compareSpectraAndStack(nameRoot):
	summaryTableName =join(stackDir, nameRoot + "summaryTable_stack.fits")
	table = fits.open(summaryTableName)[1].data
	chi1 = n.empty(len(table['PLATE']))
	chi2 = n.empty(len(table['PLATE']))
	chi3 = n.empty(len(table['PLATE']))

	for ii in range(len(table['PLATE'])):
		entry = table[ii]
		print entry
		chi1[ii], chi2[ii], chi3[ii] = compareSpectrumToStack(entry, nameRoot=nameRoot, ggrid  = ggrid, rzgrid = rzgrid, grgrid = grgrid)

	summaryTableName =join(stackDir, nameRoot + "summaryTable_stack_comparison.fits")
	col_chi1 = fits.Column(name="chi2_Z1",format="D", array= chi1)
	col_chi2 = fits.Column(name="chi2_Z2",format="D", array= chi2)
	col_chi3 = fits.Column(name="chi2_Z3",format="D", array= chi3)
	cols = table.columns + col_chi1 + col_chi2 + col_chi3
	tbhdu = fits.BinTableHDU.from_columns(cols)
	prihdr = fits.Header()
	prihdr['chunk'] = nameRoot
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	os.system('rm '+summaryTableName)
	thdulist.writeto(summaryTableName)

nameRoot="elg270_eboss17"
compareSpectraAndStack(nameRoot)
nameRoot="elg270_eboss67"
compareSpectraAndStack(nameRoot)
