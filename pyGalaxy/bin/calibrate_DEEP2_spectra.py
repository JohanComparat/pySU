"""
This script performs the flux calibration of the DEEP2 spectra
"""

from GalaxySpectrumDEEP2 import *

# defines the galaxy survey class : loads catalogs and calibration files
gs=GalaxySurveyDEEP2(redshift_catalog="zcat.deep2.dr4.fits", calibration=True)

# loops over elements of the catalog looking for the spectra and flux calibrates them

for ii in range(len(gs.catalog)):
	catalog_entry = gs.catalog[ii]
	spec1d = GalaxySpectrumDEEP2(catalog_entry,calibration=True,lineFits=False)
	if len(spec1d.path_to_spectrum)==1 and spec1d.catalog_entry['ZBEST']>=0 and spec1d.catalog_entry ['ZQUALITY']>=1 :
		print spec1d.path_to_spectrum
		spec1d.openObservedSpectrum()
		if n.max(spec1d.lambd)>6500 :
			checkList=glob.glob(spec1d.path_to_spectrum[0][:-5]+"_fc.dat")
			if len( checkList) >=1 :
				print "already FC:",spec1d.path_to_spectrum[0]
				continue

			print spec1d.slit, spec1d.mask, spec1d.ID,time.time()
			spec1d.correctQE()
			fc=spec1d.fluxCal()
			if fc=="bad":
				print "not calibrable:",spec1d.path_to_spectrum[0]
			else:
				spec1d.writeFCspec()
