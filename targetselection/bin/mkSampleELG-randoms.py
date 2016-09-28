import os
from os.path import join
import astropy.io.fits as fits
import astropy.units as u
import astropy.coordinates as co
from astropy.coordinates import SkyCoord
import numpy as n
import mangle
import sys

chunk = "eboss8"

outputDir = join(os.environ['EBOSS_ROOT'],"target","elg","clusteringSample")
outputCatalog = "ELG-catalog-angular.random.dat"

# loads geometry information :
tilingDir = join(os.environ['EBOSS_ROOT'],"ebosstilelist/trunk/outputs",chunk)
maskGeometry = mangle.Mangle(join(tilingDir,"geometry-"+chunk+".fits"))
maskGeometryTSR = mangle.Mangle(join(outputDir,"geometry-"+chunk+"-TSR.ply"))

# loads randoms that account for the geometry of the photometric survey :
raRinter, decRinter = n.loadtxt(join(os.environ['EBOSS_ROOT'],"target","elg", "decalsSelections", "randoms", "randoms_ELG_280.random.masked"), usecols=(0,1), unpack=True)
selection = (raRinter>maskGeometry.ramid.min()-3.) & (raRinter<maskGeometry.ramid.max()+3.) & (decRinter>maskGeometry.decmid.min()-3.) & (decRinter<maskGeometry.decmid.max()+3.)

raR = raRinter[selection]
decR = decRinter[selection]

pID_R = n.array([ maskGeometryTSR.polyid( raR[ii], decR[ii]) for ii in range(len(raR)) ])

raR2 ,decR2, tsrR = [], [], []
print maskGeometryTSR.weights
for index,id in enumerate(maskGeometryTSR.polyids) :
	NR = len((pID_R == id).nonzero()[0])
	rds = n.random.uniform(0,1,NR)
	sel = (rds < maskGeometryTSR.weights[index])
	raR2.append ( raR[(pID_R == id)][(sel)] )
	decR2.append ( decR[(pID_R == id)][(sel)] )
	tsrR.append ( n.ones_like(decR[(pID_R == id)][(sel)])*maskGeometryTSR.weights[index] )


# print "11"
n.savetxt(join(outputDir,outputCatalog), n.transpose([n.hstack((raR2)), n.hstack((decR2)), n.hstack((tsrR))]), header = " ra dec tsr" )

