#! /usr/bin/env python

"""
This script computes the line luminosities from of the VVDS survey
"""

# line list to be converted to luminosities
from os.path import join
import astropy.cosmology as co
cosmo=co.Planck15 #co.FlatLambdaCDM(H0=70,Om0=0.3)
import scipy.spatial.ckdtree as t

from GalaxySurveyVVDS import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import pyregion 
import pyregion._region_filter as filter
import numpy as n

# define field number 02 for the DEEP and 10, 14, 22 for the WIDE.
fields = [10, 14, 22] 
magLim =[ 22.5, 22.5, 22.5]
summaryCat = [ "VVDS_WIDE_F10_summary.LFcatalog.Planck15.v2.fits", "VVDS_WIDE_F14_summary.LFcatalog.Planck15.v2.fits", "VVDS_WIDE_F22_summary.LFcatalog.Planck15.v2.fits"]
summaryCatOut = [ "VVDS_WIDE_F10_summary.LFcatalog.Planck15.v3.fits", "VVDS_WIDE_F14_summary.LFcatalog.Planck15.v3.fits", "VVDS_WIDE_F22_summary.LFcatalog.Planck15.v3.fits"]
photozCat = 'photozCFHTLS-W1234-g25-Idef.fits'
depth = [ "WIDE", "WIDE", "WIDE"]
dBin = 0.05
bins = n.arange(0,1.4,dBin)
finalCommandConcatenate = """java -jar ~/stilts.jar tcat ifmt=fits in="VVDS_WIDE_F10_summary.LFcatalog.Planck15.v3.fits VVDS_WIDE_F14_summary.LFcatalog.Planck15.v3.fits VVDS_WIDE_F22_summary.LFcatalog.Planck15.v3.fits" out=VVDS_WIDE_summary.LFcatalog.Planck15.v3.fits"""


########################################
########################################
# SSR relation using photoz
########################################
########################################

ii = 2
#for ii in range(len(fields)):
print fields[ii]
# loads the catalog for the field
# spectroscopy + mask
survey = GalaxySurveyVVDS( redshift_catalog = summaryCat[ii] )
survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
field = (survey.catalog['NUM']/1e7).astype(int)
allSpec = (field==fields[ii])
speccat = survey.catalog[allSpec]
Nspec_total = len(speccat)

SSR = n.ones_like(speccat['ALPHA'])*-1.
SSR_ERR = n.ones_like(speccat['ALPHA'])*-1.

redshiftBest = n.empty(Nspec_total)

goodZ = (speccat['ZFLAGS']==2)|(speccat['ZFLAGS']==3)|(speccat['ZFLAGS']==4)|(speccat['ZFLAGS']==9)
badZ = (goodZ==False)

redshiftBest[goodZ] = speccat['Z'][goodZ]

raBad, decBad = speccat['ALPHA'][badZ], speccat['DELTA'][badZ]

# photometric redshift catalog
zphcatall = fits.open(join(survey.vvds_photo_dir, photozCat) )[1].data
zphcat = zphcatall[ (zphcatall['zphot_T07']>-0.1)&(zphcatall['zphot_T07']<2.) ]
treePhZ = t.cKDTree( n.transpose([zphcat['alpha_T07'], zphcat['delta_T07']]) ,1000.0)

indexes = treePhZ.query(n.transpose([raBad,decBad]),1)
dist = indexes[0]
ids = indexes[1]
ok = (dist < 1.5/3600.)
nok = (dist >= 1.5/3600.)
#redshiftBest[badZ] = zphcat['zphot_T07'][ids]
redshiftBest[badZ[ok]] = zphcat['zphot_T07'][ids[ok]]
redshiftBest[badZ[nok]] = n.ones_like(zphcat['zphot_T07'][ids[nok]])*-1.

print float(len(zphcat['zphot_T07'][ids[ok]]))/len(raBad), "% of photoz used"
print len(zphcat['zphot_T07'][ids[ok]])
print len(raBad)
print len(redshiftBest[badZ[ok]])
print len(redshiftBest[badZ[nok]])
# print ND, NR
NDpb = n.histogram(redshiftBest[goodZ], bins=bins)[0]
NRpb = n.histogram(redshiftBest, bins=bins)[0]
# piecewise interpolation of the TSR :
ssr_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), NDpb.astype(float)/NRpb)
ssr_err_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), NDpb**(-0.5) * NDpb.astype(float)/NRpb)

SSR = ssr_eval(redshiftBest[goodZ])
SSR_ERR = ssr_err_eval(redshiftBest[goodZ])

# writes the new catalog
c0 = fits.Column(name="SSR_new",format="D", array= SSR )
c1 = fits.Column(name="SSR_ERR_new",format="D", array= SSR_ERR )
new_columns = survey.catalog.columns + c0 + c1
hdu = fits.BinTableHDU.from_columns(new_columns)
os.system("rm -rf "+join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]) )
hdu.writeto(join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]))


########################################
########################################
# extrapolation to other field
########################################
########################################

ii = 0
#for ii in range(len(fields)):
print fields[ii]
# loads the catalog for the field
# spectroscopy + mask
survey = GalaxySurveyVVDS( redshift_catalog = summaryCat[ii] )
survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
field = (survey.catalog['NUM']/1e7).astype(int)
allSpec = (field==fields[ii])
speccat = survey.catalog[allSpec]
Nspec_total = len(speccat)

SSR = n.ones_like(speccat['ALPHA'])*-1.
SSR_ERR = n.ones_like(speccat['ALPHA'])*-1.

redshiftBest = speccat['Z']

goodZ = (speccat['ZFLAGS']==2)|(speccat['ZFLAGS']==3)|(speccat['ZFLAGS']==4)|(speccat['ZFLAGS']==9)
badZ = (goodZ==False)

SSR = ssr_eval(redshiftBest[goodZ])
SSR_ERR = ssr_err_eval(redshiftBest[goodZ])


c0 = fits.Column(name="SSR_new",format="D", array= SSR )
c1 = fits.Column(name="SSR_ERR_new",format="D", array= SSR_ERR )
new_columns = survey.catalog.columns + c0 + c1
hdu = fits.BinTableHDU.from_columns(new_columns)
os.system("rm -rf "+join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]) )
hdu.writeto(join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]))


ii = 1
#for ii in range(len(fields)):
print fields[ii]
# loads the catalog for the field
# spectroscopy + mask
survey = GalaxySurveyVVDS( redshift_catalog = summaryCat[ii] )
survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
field = (survey.catalog['NUM']/1e7).astype(int)
allSpec = (field==fields[ii])
speccat = survey.catalog[allSpec]
Nspec_total = len(speccat)

SSR = n.ones_like(speccat['ALPHA'])*-1.
SSR_ERR = n.ones_like(speccat['ALPHA'])*-1.

redshiftBest = speccat['Z']

goodZ = (speccat['ZFLAGS']==2)|(speccat['ZFLAGS']==3)|(speccat['ZFLAGS']==4)|(speccat['ZFLAGS']==9)
badZ = (goodZ==False)

SSR = ssr_eval(redshiftBest[goodZ])
SSR_ERR = ssr_err_eval(redshiftBest[goodZ])


c0 = fits.Column(name="SSR_new",format="D", array= SSR )
c1 = fits.Column(name="SSR_ERR_new",format="D", array= SSR_ERR )
new_columns = survey.catalog.columns + c0 + c1
hdu = fits.BinTableHDU.from_columns(new_columns)
os.system("rm -rf "+join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]) )
hdu.writeto(join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]))

os.system( finalCommandConcatenate )

"""
Eventually toimplement per small units of mask, but beware of overlaps between masks

tsr_A = n.empty((len(spec_mask_filter), len(bins)-1))
nInReg = n.empty((len(spec_mask_filter), 3 ))
AREA = n.zeros_like(speccat['ALPHA'])
TSR = n.zeros_like(speccat['ALPHA'])
TSR_ERR = n.zeros_like(speccat['ALPHA'])
randRA = []
randDEC = []
for jj, mask in enumerate(spec_mask_filter):
	areaIn = mask.inside(n.transpose([raPix, decPix]))
	NPIX = len(areaIn.nonzero()[0])
	specIn = mask.inside(n.transpose([speccat['ALPHA'], speccat['DELTA']]))
	ND = len(specIn.nonzero()[0])
	photIn = mask.inside(n.transpose([photocat['ALPHA'], photocat['DELTA']]))
	NR = len(photIn.nonzero()[0])
	randIn = mask.inside(n.transpose([randRA_all, randDEC_all]))
	downSelectRD = (n.random.uniform(0,1,len(randRA_all[randIn]))< float(ND)/float(NR) )
	randRA.append( randRA_all[randIn][downSelectRD] )
	randDEC.append( randDEC_all[randIn][downSelectRD] )
	# print ND, NR
	NDpb = n.histogram(speccat['MAGI'][specIn], bins=bins)[0]
	NRpb = n.histogram(photocat['MAGI_CFH12K'][photIn], bins=bins)[0]
	tsr_A[jj] = NDpb.astype(float)/NRpb
	tsr_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), tsr_A[jj])
	tsr_err_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), NDpb**(-0.5) *tsr_A[jj])
	TSR[specIn] = tsr_eval(speccat['MAGI'][specIn])
	TSR_ERR[specIn] = tsr_err_eval(speccat['MAGI'][specIn])
	area = NPIX * areaPerPixel  #spec_mask_polygon[jj]
	AREA[specIn] = n.ones_like(AREA[specIn]) * area
	nInReg[jj] = ND, NR, area
	# print tsr_A[jj]

randRA = n.hstack(( randRA ))
randDEC = n.hstack(( randDEC ))
areaEffective = nInReg.T[2].sum()
	
for kk in range(len(tsr_A)):
	p.plot((bins[1:]+bins[:-1])/2., tsr_A[kk])

p.savefig(join(os.environ['DATA_DIR'],'indev','tsr'+str(fields[ii])+".png"))

p.plot(photocat['ALPHA'][photIn], photocat['DELTA'][photIn], 'b,')
p.plot(speccat['ALPHA'][specIn], speccat['DELTA'][specIn], 'r,')
p.savefig('ra-dec-2'+str(fields[ii])+".png")

p.plot(photoDat['ALPHA'], photoDat['DELTA'], 'b,',label='phot')
p.plot(survey.catalog['ALPHA'][allSpec], survey.catalog['DELTA'][allSpec], 'r+', label='spec')
p.legend()
p.savefig('ra-dec'+str(fields[ii])+".png")


##############################################3
##############################################3
# writes the resutls
##############################################3
##############################################3
##############################################3exit

new_columns = survey.catalog.columns

for line in lineList :
	print line[2]
	c0,c1 = survey.computeLineLuminosity(line,sphereCM)
	new_columns += c0 
	new_columns += c1

hdu = fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)

"""