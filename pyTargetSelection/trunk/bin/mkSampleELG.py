import os
from os.path import join
import astropy.io.fits as fits
import astropy.units as u
import astropy.coordinates as co
from astropy.coordinates import SkyCoord
import numpy as n
import mangle
import sys

# parameters :
chunk = "eboss8"
dMax = 0.05 # maximum magnitude separation for a collided galaxy to be overweighted for in its nearest neighbor.

outputDir = join(os.environ['EBOSS_ROOT'],"target","elg","clusteringSample")
outputCatalog = "ELG-catalog-angular.fits"

# first, compute angular weights : TSR and fiber collisions

# opens the photometry from which the targets were issued
print "computes the target sampling rate in eboss8 ELG chunk, takes into account fiber collisions"
hduPhotoObjAll = fits.open(join(os.environ['EBOSS_ROOT'],"target","elg","decalsSelections","catalog_ELG_280-masked-readyForTiling.fits"))

brickid = hduPhotoObjAll[1].data['brickid']
objid = hduPhotoObjAll[1].data['objid']
strID = n.array([ str(brickid[ii])+"_"+str(objid[ii]) for ii in range(len(objid)) ])
ra = hduPhotoObjAll[1].data['ra']
dec = hduPhotoObjAll[1].data['dec']

g = 22.5 - 2.5 * n.log10(hduPhotoObjAll[1].data['decam_flux'].T[1] / hduPhotoObjAll[1].data['decam_mw_transmission'].T[1])
r = 22.5 - 2.5 * n.log10(hduPhotoObjAll[1].data['decam_flux'].T[2] / hduPhotoObjAll[1].data['decam_mw_transmission'].T[2])
z = 22.5 - 2.5 * n.log10(hduPhotoObjAll[1].data['decam_flux'].T[4] / hduPhotoObjAll[1].data['decam_mw_transmission'].T[4])

# tiling information
# print "2"
tilingDir = join(os.environ['EBOSS_ROOT'],"ebosstilelist/trunk/outputs",chunk)
hduTiled = fits.open(join(tilingDir,"final-"+chunk+".fits"))
#hduTiled = fits.open(join(tilingDir,"all-"+chunk+".fits"))
raTiled = hduTiled[1].data['RA']
decTiled = hduTiled[1].data['DEC']

brickidTiled = hduTiled[1].data['BRICKID_DECALS']
objidTiled = hduTiled[1].data['OBJID_DECALS']
strIDTiled = n.array([ str(brickidTiled[ii])+"_"+str(objidTiled[ii]) for ii in range(len(objidTiled)) ])

# print "strIDTiled", strIDTiled, len(strIDTiled)

maskGeometry = mangle.Mangle(join(tilingDir,"geometry-"+chunk+".ply"))

observed = n.zeros_like(strID,dtype=bool)
for jj in range(len(strID)):
	if len((strID[jj] == strIDTiled).nonzero()[0])>0 :
		observed[jj] = True 
	else :
		observed[jj] =False

len(observed.nonzero()[0])
# print "3"
# compute TSR : N galaxy with fiber / N galaxy targeted in each polygon.
pID = n.array([ maskGeometry.polyid( ra[ii], dec[ii]) for ii in range(len(ra)) ])
pID_tiled = n.array([ maskGeometry.polyid( raTiled[ii], decTiled[ii]) for ii in range(len(raTiled)) ])

fcRadius = 62. # arcseconds
TSR = n.empty(len(raTiled))
weight = n.ones_like(raTiled) # ones everywhere

for index,id in enumerate(maskGeometry.polyids) :
	Ntargeted = len((pID == id).nonzero()[0])
	NwithFiber = len((pID_tiled == id).nonzero()[0])
	print NwithFiber, Ntargeted 
	if Ntargeted == NwithFiber : # exactly complete
		TSR[(pID_tiled == id)] = n.ones_like(TSR[(pID_tiled == id)])
		maskGeometry.weights[index] = 1.
		# print "5"
	else : # not complete, possible fiber collisions
		# defines the targets and fiber assigned targets in the polygon and that did not get a fiber
		selT = (pID == id)&(observed==False)
		raT,decT=ra[selT],dec[selT]
		selF = (pID_tiled == id)
		raF,decF,strIDTiledF=raTiled[selF],decTiled[selF],strIDTiled[selF]
		# print raF,decF,strIDTiledF,len(raF),len(decF),len(strIDTiledF)
		cT = SkyCoord(raT, decT, frame='icrs', unit='deg')
		cF = SkyCoord(raF, decF, frame='icrs', unit='deg')
		# print "6"
		# loops over the observed objects if a non-observed target was fiber collided
		for k in range(len(raT)):
			collided = (cF.separation(cT[k]) < fcRadius *u.deg/3600.)
			Ncollided = len((cF.separation(cT[k]) < fcRadius *u.deg/3600.).nonzero()[0])
			# print collided, Ncollided
			# print "61"
			if 	Ncollided == 0 : # not complete, no fiber collision
				TSR[(pID_tiled == id)] = n.ones_like(TSR[(pID_tiled == id)])* float(NwithFiber) /float(Ntargeted)
				maskGeometry.weights[index] = float(NwithFiber) /float(Ntargeted)
				# print "62"
			if 	Ncollided > 0 : # not complete, there is fiber collision ...
				# print "Ncollided", Ncollided
				# print strIDTiledF
				count = 0
				for i in range(Ncollided):
					# print strIDTiledF[collided][i]
					# print len(strIDTiledF), len(collided), len(weight), len(strIDTiled)
					# photometry for both galaxy :
					distanceMagnitude = ( abs(g[(pID == id)][(collided)][i] - g[(strID == strIDTiledF[collided][i])])< dMax) & ( abs(r[(pID == id)][(collided)][i] - r[(strID == strIDTiledF[collided][i])])< dMax ) & ( abs(z[(pID == id)][(collided)][i] - z[(strID == strIDTiledF[collided][i])])<dMax ) 
					# print distanceMagnitude
					if distanceMagnitude[0] :
						weight[(strIDTiled == strIDTiledF[collided][i])] = 2
						#weight[selF][collided][i] = 2
						count += 1

				# print "7"
				TSR[(pID_tiled == id)] = n.ones_like(TSR[(pID_tiled == id)])* float(NwithFiber + count) /float(Ntargeted)				
				maskGeometry.weights[index] = float(NwithFiber+count) / float(Ntargeted)

# We write the result in a catalog
# print "8"
tsrCol = fits.Column(name="TSR",format="D", array=TSR )
MaskCol = fits.Column(name="weight",format="D", array=weight )
newcols = hduTiled[1].data.columns + tsrCol + MaskCol

hdu = fits.BinTableHDU.from_columns(newcols)
os.system("rm -rf "+join(outputDir,outputCatalog))
hdu.writeto(join(outputDir,outputCatalog))

maskGeometry.writeply(join(outputDir,"geometry-"+chunk+"-TSR.ply"))


