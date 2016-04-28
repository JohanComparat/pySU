#! /usr/bin/env python

"""
This script computes the line luminosities from of the VVDS survey
"""

# line list to be converted to luminosities
from os.path import join
import astropy.cosmology as co
cosmo=co.Planck15 #co.FlatLambdaCDM(H0=70,Om0=0.3)

from GalaxySurveyVVDS import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import pyregion 
import pyregion._region_filter as filter
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
	return mask_filter

myfilter = filter.Polygon(36.5232083,-4.5198276,36.4187083,-4.5198276,36.4089233,-4.5270616,36.4089233,-4.6536946,36.5232083,-4.6536946,36.5232083,-4.5198276) | filter.Polygon(36.3787413,-4.5289276,36.3749413,-4.5219276,36.2662083,-4.5219276,36.2662083,-4.6514946,36.3787413,-4.6514946,36.3787413,-4.5289276)
myfilter = filter.Polygon(878.11237,811.88766,1003,861,973.00003,688,878.11237,752.11234)

myfilter = filter.Polygon([43.1, 43.5, 43.6],[0.5, 0.8,0.5])

myfilter.inside1(43.3, 0.501)
myfilter.inside1(10, 0)
myfilter.inside([0, 10], [0, 0])

# define field number 02 for the DEEP and 10, 14, 22 for the WIDE.
fields = [02, 10, 14, 22] 
magLim =[24., 22.5, 22.5, 22.5]
summaryCat = ["VVDS_DEEP_summary.LFcatalog.Planck15.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.fits"]
summaryCatOut = ["VVDS_DEEP_summary.LFcatalog.Planck15.v2.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.v2.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.v2.fits", "VVDS_WIDE_summary.LFcatalog.Planck15.v2.fits"]
depth = ["DEEP", "WIDE", "WIDE", "WIDE"]
dMag = 0.5

ii = 0

# loads the catalog for the field
# spectroscopy + mask
survey = GalaxySurveyVVDS( redshift_catalog = summaryCat[ii] )
survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
field = (survey.catalog['NUM']/1e7).astype(int)
allSpec = (field==fields[ii])
speccat = survey.catalog[allSpec]
path_to_spec_mask = join(survey.vvds_catalog_dir, "F"+str(fields[ii]).zfill(2)+"_specmask.reg")

# photometry
# catalog of galaxies + mask with the VVDS selection applied
path_to_photo_Cat = join(survey.vvds_photo_dir, "cesam_vvds_photoF"+str(fields[ii]).zfill(2)+"_"+depth[ii]+".fits")
photocat_all = fits.open(path_to_photo_Cat)[1].data
selection=(photocat_all['MAGI']<magLim[ii])&(photocat_all['MAGI']>17.5)
photocat = photocat_all[selection]
path_to_photo_mask = join(survey.vvds_photo_dir, "F"+str(fields[ii]).zfill(2)+"_photmask.reg")

# loads the actual co-added photo NOT USED SO FAR
path_to_photo_Iband = join(survey.vvds_photo_dir, "F"+str(fields[ii]).zfill(2)+"i.fits")
fitsfile = fits.open(path_to_photo_Iband)
iband_header = fitsfile[0].header
iband_data = fitsfile[0].data

# creates the masks
spec_mask_filter = createMask(path_to_spec_mask)
photo_mask_filter = createMask(path_to_photo_mask)

bins = n.arange(17, magLim[ii], dMag)

tsr_A = n.empty((len(spec_mask_filter), len(bins)-1))
for jj, mask in enumerate(spec_mask_filter):
	specIn = mask.inside(n.transpose([speccat['ALPHA'], speccat['DELTA']]))
	ND = len(specIn.nonzero()[0])
	photIn = mask.inside(n.transpose([photocat['ALPHA'], photocat['DELTA']]))
	NR = len(specIn.nonzero()[0])
	NDpb = n.histogram(speccat['MAGI'][specIn], bins=bins)[0]
	NRpb = n.histogram(photocat['MAGI'][photIn], bins=bins)[0]
	tsr_A[jj] = float(NDpb)/float(NRpb)
	print ND, NR, tsr_A[jj]

"""
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