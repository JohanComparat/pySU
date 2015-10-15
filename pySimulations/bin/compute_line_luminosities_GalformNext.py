"""
This script computes the line luminosities from of the Galform survey
"""
# line list to be converted to luminosities
from lineListAir import *
lineList = n.array([ [O2,O2_mean,"O2_3728"], [O3,O3_5007,"O3_5007"], [H1,H1_4862,"H1_4862"] ])

import astropy.cosmology as co
cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)

from GalaxySurveyGalform import *

# Galform catalog
survey = GalaxySurveyGalform(redshift_catalog = "galform.elg.next.fits")
survey.path_to_output_catalog  =  survey.galform_catalog_dir+ "galform.next.LFcatalog.fits"
survey.Ngalaxies=len(survey.catalog)
new_columns = survey.catalog.columns

dL=cosmo.luminosity_distance(survey.catalog['zObs'])
sphere=4*n.pi*(dL)**2.
sphereCM=sphere.to(u.cm**2)

for line in lineList :
	print line[2]
	c0 = survey.computeLineLuminosity(line,sphereCM)
	new_columns += c0 

hdu = fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)
