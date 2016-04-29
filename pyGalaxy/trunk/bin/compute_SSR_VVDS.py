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
import healpy as hp
NSIDE = 4096
ids = n.arange( hp.nside2npix( NSIDE ) )
areaPerPixel = 129600. / n.pi / len( ids )
decPixrad, raPixrad = hp.pix2ang(NSIDE, ids )
raPixAll = raPixrad*180./n.pi
decPixAll = decPixrad*180./n.pi -90.

# creates the filter masks
def createMask(path_to_mask):
	"""
	From a region file containing only polygons, it creates a maskfile with "inside" method to filter the data.
	"""
	f=open(path_to_mask,'r')
	lines = n.array(f.readlines())
	f.close()
	polygons = []
	for ii,line in enumerate(lines):
		if line[:3]=="pol":
			arr = n.array(line[8:-2].split(',')).astype('float')
			polygons.append(arr.reshape(len(arr)/2,2))

	filters = [ (filter.Polygon(polygons[0].T[0], polygons[0].T[1])) | (filter.Polygon(polygons[1].T[0], polygons[1].T[1])) ]
	for ii in n.arange(2,len(polygons),1) :
		filters  = [ (filters[0]) | (filter.Polygon(polygons[ii].T[0], polygons[ii].T[1])) ]

	mask_filter = filters[0]
	return mask_filter, polygons


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



ii = 0
for ii in range(len(fields)):
	print fields[ii]
	# loads the catalog for the field
	# spectroscopy + mask
	survey = GalaxySurveyVVDS( redshift_catalog = summaryCat[ii] )
	survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
	field = (survey.catalog['NUM']/1e7).astype(int)
	allSpec = (field==fields[ii])
	speccat = survey.catalog[allSpec]
	path_to_spec_mask = join(survey.vvds_catalog_dir, "F"+str(fields[ii]).zfill(2)+"_specmask.reg")
	Nspec_total = len(speccat)

	SSR = n.zeros_like(speccat['ALPHA'])
	SSR_ERR = n.zeros_like(speccat['ALPHA'])

	redshiftBest = n.empty(Nspec_total)
	
	goodZ = (speccat['ZFLAGS']==2)|(speccat['ZFLAGS']==3)|(speccat['ZFLAGS']==4)|(speccat['ZFLAGS']==9)
	badZ = (goodZ==False)
	
	redshiftBest[goodZ] = speccat['Z'][goodZ]
	
	raBad, decBad = speccat['ALPHA'][badZ], speccat['DELTA'][badZ]
	
	# photometric redshift catalog
	zphcatall = fits.open(join(survey.vvds_photo_dir, photozCat) )[1].data
	zphcat = zphcatall [][]
	treePhZ = t.cKDTree(n.transpose([zphcat['alpha'], zphcat['delta']))#,1000.0)
	indexes = treePhZ.query([raBad,decBad],1)
	redshiftBest[badZ] = zphcat['zphot'][indexes]
	
	# print ND, NR
	NDpb = n.histogram(redshiftBest[goodZ], bins=bins)[0]
	NRpb = n.histogram(redshiftBest, bins=bins)[0]
	# piecewise interpolation of the TSR :
	tsr_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), NDpb.astype(float)/NRpb)
	tsr_err_eval = lambda x : n.piecewise(x, n.array([ (x > bins[kk])&(x<bins[kk+1]) for kk in range(len(NDpb)) ]), NDpb**(-0.5) * NDpb.astype(float)/NRpb)

	SSR = tsr_eval(redshiftBest)
	SSR_ERR = tsr_err_eval(redshiftBest)

	# writes the new catalog
	c0 = fits.Column(name="TSR_new",format="D", array= TSR )
	c1 = fits.Column(name="TSR_ERR_new",format="D", array= TSR_ERR )
	new_columns = survey.catalog.columns + c0 + c1
	hdu = fits.BinTableHDU.from_columns(new_columns)
	os.system("rm -rf "+join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]) )
	hdu.writeto(join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii]))

	# writes the randoms 
	col0 = fits.Column(name="RA",format='D', array= randRA_all[randIn]) #[downSelectRD])
	col1 = fits.Column(name="DEC",format='D', array= randDEC_all[randIn]) #[downSelectRD])
	fullSpec_cols  = fits.ColDefs([col0, col1])
	fullSpec_tb_hdu = fits.BinTableHDU.from_columns(fullSpec_cols)
	os.system("rm -rf "+join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii][:-5]+".random.fits") )
	fullSpec_tb_hdu.writeto(join(os.environ['VVDS_DIR'], 'catalogs', summaryCatOut[ii][:-5]+".random.fits") )

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