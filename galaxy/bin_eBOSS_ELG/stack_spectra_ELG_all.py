#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"SDSS/stacks")
specList = join(spec_dir, "eboss-elg_0.2_z_1.5.asc") 
outfile = join(spec_dir, os.path.basename(specList)[:-4]+".specMatrix")
stack=sse.SpectraStackingEBOSS(specList, outfile)

def getSpectra(path_to_spectrum):
	hdulist = fits.open(path_to_spectrum)
	wave = 10**hdulist[1].data['loglam']
	ok=(wave>3740)&(wave<9604)
	flux = hdulist[1].data['flux']
	hdulist.close()
	return wave[ok], flux[ok]

for IDX_j in n.arange(0, len(stack.plates), 4096):
	IDX_min=IDX_j
	IDX_max=IDX_j+4096
	IDX_str=str(IDX_min).zfill(6)+'-'+str(IDX_max).zfill(6)
	samp_plates, samp_mjds, samp_fiberids, samp_redshifts = stack.plates[IDX_min:IDX_max], stack.mjds[IDX_min:IDX_max], stack.fiberids[IDX_min:IDX_max], stack.redshifts[IDX_min:IDX_max]

	FLUXES = n.zeros((samp_plates.shape[0], 4096))
	data = []
	bad_ids = []
	for jj, (plate, mjd, fiber, redshift) in enumerate(zip( samp_plates, samp_mjds, samp_fiberids, samp_redshifts )):
		path_to_spectrum = sse.get_path_to_spectrum_v5_11_0(plate, mjd, fiber)
		if os.path.isfile(path_to_spectrum):
			wl,fl=getSpectra(path_to_spectrum)
			data.append([fl.shape[0], wl.min(), wl.max()])
			if fl.shape[0]==4096:
				FLUXES[jj]=fl
				wavelength=wl

	n.savetxt(stack.out_file+'.'+IDX_str+'.dat', FLUXES)
	n.savetxt(stack.out_file+'.wavelength.'+IDX_str+'.dat', wavelength)
	n.savetxt(stack.out_file+'.shapes.'+IDX_str+'.dat', n.array(data) )
	n.savetxt(stack.out_file+'.list.'+IDX_str+'.dat', n.array([samp_plates, samp_mjds, samp_fiberids, samp_redshifts]) )

