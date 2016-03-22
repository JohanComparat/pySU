import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
from HandleSdssPlate import *
import glob	
specDirSdssMain = "/uufs/chpc.utah.edu/common/home/sdss00/sdsswork/sdss/spectro/redux/26/spectra"
list_of_stacks_eb67 = glob.glob(join("/uufs/chpc.utah.edu/common/home/u0992342/eboss67/grz_stacks/","*.asc"))
list_of_stacks_eb17 = glob.glob(join("/uufs/chpc.utah.edu/common/home/u0992342/eboss17/grz_stacks/","*.asc"))

list_of_stacks = n.hstack((list_of_stacks_eb67,list_of_stacks_eb17))
list_of_stacks.sort()

for ii, el in enumerate(list_of_stacks):
	PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol = n.loadtxt(el, unpack=True, usecols=(0,1,2,3,4,5,6))
	g_min = n.min(gmag)
	g_max = n.max(gmag)
	gr_min = n.min(grcol)
	gr_max = n.max(grcol)
	rz_min = n.min(grcol)
	rz_max = n.max(grcol)
	stackName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG",el[:-4].split('/')[-1] + "_stack.fits")
	outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG",el[:-4].split('/')[-1] + "_stackComparison.dat")
	hdu = fits.open(stackName)[1].data
	sel= (hdu['NspectraPerPixel']>0.9*n.max(hdu['NspectraPerPixel']))
	chi2median = n.empty_like(PLATE) 
	chi2mean = n.empty_like(PLATE) 
	for jj in range(len(PLATE)):
		ObsPlate = HandleReducedELGPlate(int(PLATE[ii]),int(MJD[ii]))
		ObsPlate.loadSpec(int(FIBERID[ii]))
		wlmin=n.min(hdu['wavelength'][sel]*(1+REDSHIFT[jj]))
		wlmax=n.max(hdu['wavelength'][sel]*(1+REDSHIFT[jj]))
		medianStack = interp1d(hdu['wavelength'][sel]*(1+REDSHIFT[jj]),hdu['meanWeightedStack'][sel])
		meanStack =interp1d(hdu['wavelength'][sel]*(1+REDSHIFT[jj]),hdu['medianStack'][sel])
		overlap=(ObsPlate.wavelength>wlmin)&(ObsPlate.wavelength<wlmax)
		x = ObsPlate.wavelength[overlap]
		y = ObsPlate.flux[overlap]
		yerr = ObsPlate.fluxErr[overlap]
		chi2median[jj] = n.sum(((y - medianStack(x))/yerr)**2.)/len(x)
		chi2mean[jj] = n.sum(((y - meanStack(x))/yerr)**2.)/len(x)

	header = " PLATE MJD FIBERID REDSHIFT gmag rzcol grcol chi2median chi2mean"
	n.savetxt(outPutFileName, n.transpose([PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol , chi2median, chi2mean]), header = header)

