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

myfilter = filter.Polygon(36.5232083,-4.5198276,36.4187083,-4.5198276,36.4089233,-4.5270616,36.4089233,-4.6536946,36.5232083,-4.6536946,36.5232083,-4.5198276) | filter.Polygon(36.3787413,-4.5289276,36.3749413,-4.5219276,36.2662083,-4.5219276,36.2662083,-4.6514946,36.3787413,-4.6514946,36.3787413,-4.5289276)
myfilter = filter.Polygon(878.11237,811.88766,1003,861,973.00003,688,878.11237,752.11234)

myfilter = filter.Polygon([43.1, 43.5, 43.6],[0.5, 0.8,0.5])

myfilter.inside1(43.3, 0.501)
myfilter.inside1(10, 0)
myfilter.inside([0, 10], [0, 0])

# define field number 02 for the DEEP and 10, 14, 22 for the WIDE.
fields = [02, 10, 14, 22] 
ii = 0

# loads the catalog for the VVDS-Deep
# spectroscopy
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_DEEP_summary.LFcatalog.Planck15.fits")
survey.vvds_photo_dir = join(os.environ['VVDS_DIR'],'photometry')
field = (survey.catalog['NUM']/1e7).astype(int)
galaxies= n.transpose([survey.catalog['ALPHA'],survey.catalog['DELTA']])
allSpec = (field==fields[ii])
path_to_spec_mask = join(survey.vvds_catalog_dir, "F"+str(fields[ii]).zfill(2)+"_specmask.reg")
speccat = survey.catalog

# photometry
# catalog of galaxies
path_to_photo_Cat = join(survey.vvds_photo_dir, "cesam_vvds_photoF"+str(fields[ii]).zfill(2)+"_DEEP.fits")
photocat = fits.open(path_to_photo_Cat)[1].data
# actual co-added photo
path_to_photo_Iband = join(survey.vvds_photo_dir, "F"+str(fields[ii]).zfill(2)+"i.fits")
fitsfile = fits.open(path_to_photo_Iband)
iband_header = fitsfile[0].header
iband_data = fitsfile[0].data
path_to_photo_mask = join(survey.vvds_photo_dir, "F"+str(fields[ii]).zfill(2)+"_photmask.reg")

# creates the masks
f=open(path_to_spec_mask,'r')
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

spec_mask_filter = filters[0]
specIn = spec_mask_filter.inside(n.transpose([speccat['ALPHA'], speccat['DELTA']]))
photIn = spec_mask_filter.inside(n.transpose([photocat['ALPHA'], photocat['DELTA']]))

p.plot(photocat['ALPHA'][photIn], photocat['DELTA'][photIn], 'b,')
p.plot(speccat['ALPHA'][specIn], speccat['DELTA'][specIn], 'r,')
p.savefig('ra-dec-2'+str(fields[ii])+".png")

p.plot(photoDat['ALPHA'], photoDat['DELTA'], 'b,',label='phot')
p.plot(survey.catalog['ALPHA'][allSpec], survey.catalog['DELTA'][allSpec], 'r+', label='spec')
p.legend()
p.savefig('ra-dec'+str(fields[ii])+".png")


##############################################3
##############################################3
# convert filtersto a mask 
##############################################3
##############################################3
##############################################3exit
	
prm = n.array(rspec[0].params)

myfilterspec = rspec.get_filter()
myfilterspec.inside(galaxies)
myfilterspec.inside1(galaxies[0][0],galaxies[0][1])

rphot = pyregion.open(path_to_photo_mask).as_imagecoord(header=iband_header)
myfilterphot = rphot.get_filter()
inc = myfilterphot.inside(galaxies)
myfilterphot.inside1(galaxies[0][0],galaxies[0][1])


regL = pyregion.ShapeList([])
for ii,reg in enumerate(reglist[:10]):
	print reg
	regL.append(reg)

r2 = pyregion.parse(reglist).as_imagecoord(fitsfile[0].header)


regmask = regL.get_mask(hdu=fitsfile[0])
  


r = reg.open(path_to_photo_mask)

myfilter = filter.Polygon(r[jj].params[0], r[jj].params[1], r[jj].params[2], r[jj].params[3])

myfilter = filter.Polygon((150., 1.9), (151., 1.9), (150., 2.1)) 

(filter.Polygon =r[jj])& (filter.Polygon =r[jj+1])


jj = 0
ff = filter.r[jj]

filter = r[0] | r[1]

import pyregion._region_filter as filter
>>> myfilter = filter.Circle(0, 0, 10) & filter.Box(15, 0, 10, 10)
>>> myfilter.inside1(0, 0)
0
>>> myfilter.inside1(10, 0)
1
>>> myfilter.inside([0, 10], [0, 0])

photoDat = fits.open(join(survey.vvds_photo_dir, "cesam_vvds_photoF"+str(fields[ii])+"_WIDE.fits" ))[1].data

len(survey.catalog['NUM'][allSpec])
len(photoDat['NUM'])

p.plot(photoDat['ALPHA'], photoDat['DELTA'], 'b,',label='phot')
p.plot(survey.catalog['ALPHA'][allSpec], survey.catalog['DELTA'][allSpec], 'r+', label='spec')
p.legend()
p.savefig('ra-dec'+str(fields[ii])+".png")

survey.path_to_output_catalog  =  join(survey.vvds_catalog_dir, "VVDS_WIDE_summary.LFcatalog.Planck15.fits")
survey.Ngalaxies=len(survey.catalog)
new_columns = survey.catalog.columns

dL=cosmo.luminosity_distance(survey.catalog['Z'])
sphere=4*n.pi*(dL)**2.
sphereCM=sphere.to(u.cm**2)

for line in lineList :
	print line[2]
	c0,c1 = survey.computeLineLuminosity(line,sphereCM)
	new_columns += c0 
	new_columns += c1

hdu = fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)

# VVDS UDEEP catalog
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_UDEEP_summary.linesFitted.fits")
survey.path_to_output_catalog  =  join(survey.vvds_catalog_dir, "VVDS_UDEEP_summary.LFcatalog.Planck15.fits")
survey.Ngalaxies=len(survey.catalog)
new_columns = survey.catalog.columns

dL=cosmo.luminosity_distance(survey.catalog['Z'])
sphere=4*n.pi*(dL)**2.
sphereCM=sphere.to(u.cm**2)

for line in lineList :
	print line[2]
	c0,c1 = survey.computeLineLuminosity(line,sphereCM)
	new_columns += c0 
	new_columns += c1

hdu = fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)

# VVDS DEEP catalog
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_DEEP_summary.linesFitted.fits")
survey.path_to_output_catalog  =  join(survey.vvds_catalog_dir, "VVDS_DEEP_summary.LFcatalog.Planck15.fits")
survey.Ngalaxies=len(survey.catalog)
new_columns = survey.catalog.columns

dL=cosmo.luminosity_distance(survey.catalog['Z'])
sphere=4*n.pi*(dL)**2.
sphereCM=sphere.to(u.cm**2)

for line in lineList :
	print line[2]
	c0,c1 = survey.computeLineLuminosity(line,sphereCM)
	new_columns += c0 
	new_columns += c1

hdu = fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)
