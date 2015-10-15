"""
This script uses the LineLuminosityFunction class to estimate the LF in the emission lines for the DEEP2, VVDS and VIPERS survey.
"""

from LineLuminosityFunction import *

from lineListAir import *
import glob



# from VVDS DEEP survey
print "VVDS DEEP"
zsVIMOSmin=n.array([0.18,0.41,0.51,0.56,0.65,0.84, 1.1])
zsVIMOSmax=n.array([0.41,0.65,0.7,0.83,0.84,1.1, 1.3])
linesFittedVIMOS=n.array([ [[H1_4862,"H1_4862"],[O3_5007,"O3_5007"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O2_mean,"O2_3728"]],[[O2_mean,"O2_3728"]] ])

areaDeep=0.61
area = areaDeep
#areaUDeep=512./3600.

for ii in range(len(zsVIMOSmin)):
	zmin = zsVIMOSmin[ii]
	zmax = zsVIMOSmax[ii]
	lineSet=linesFittedVIMOS[ii]
	for line in lineSet :		
		lf = LineLuminosityFunction(lineWavelength=line[0], lineName=line[1], cosmology = cosmo, surveyName =  "VVDS", redshift_catalog = "VVDS_DEEP_summary.LFcatalog.fits", luminosityBins = n.logspace(38,45,25), Nstack = 400, Nclustering = 400, outputFolder="emissionLineLuminosityFunctions/" , zmin = zmin, zmax = zmax)
		lf.setRedshiftArray( redshiftColumn='Z' )
		lf.setRedshiftSelection( redshiftQualityColumn='ZFLAGS', lowerBound=1.9, upperBound=9.1)
		lf.setWeightArray( 1./(area * lf.catalog['SSR']*lf.catalog['TSR']), area )
		selection = (lf.catalog['TSR']>0) & (lf.catalog['SSR']>0)
		lf.computeHistogramLF(selection)
		print "---------------------------------------------------"
		print line, zmin, zmax, lf.ngals
		lf.computeHistogramVariance(selection,jk=0.1)
		lf.computeMeanWeightedRedshift(selection)
		lf.get_completness_limit(selection)
		lll=n.array( glob.glob( "/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/" + line[1] + "/" + line[1] + "-GALFORM-VVDSNEXT-z*Rcorrection.txt" ))
		redshifts=[]
		for el in lll :
			base = el.split('/')[-1]
			out = base.split('-')
			redshifts.append(float(out[3][1:5]))

		redshifts = n.array(redshifts)
		#print redshifts, lf.meanRedshift
		sel = (lf.meanRedshift - 0.05 < redshifts) & (redshifts < lf.meanRedshift + 0.05)
		#print lll[sel]
		#mCorr, error  = n.loadtxt(lll[sel][0],unpack=True ) 
		lf.correctVolumeEffect(n.ones_like(lf.LF), n.zeros_like(lf.LF))#mCorr, error)
		lf.writeLF(selection,surveyNameSuffix="DEEP")


import sys
sys.exit()

# from VVDS WIDE survey
areaWide=4.0+2.2+1.9
area = areaWide
zsVIMOSmin=n.array([0.18,0.41,0.51,0.56,0.65,0.84, 1.1])
zsVIMOSmax=n.array([0.41,0.65,0.7,0.83,0.84,1.1, 1.3])
linesFittedVIMOS=n.array([ [[H1_4862,"H1_4862"],[O3_5007,"O3_5007"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O2_mean,"O2_3728"]],[[O2_mean,"O2_3728"]] ])


for ii in range(len(zsVIMOSmin)):
	zmin = zsVIMOSmin[ii]
	zmax = zsVIMOSmax[ii]
	lineSet=linesFittedVIMOS[ii]
	for line in lineSet :
		print "---------------------------------------------------"
		print line, zmin, zmax
		lf = LineLuminosityFunction(lineWavelength=line[0], lineName=line[1], cosmology = cosmo, surveyName =  "VVDS", redshift_catalog = "VVDS_WIDE_summary.LFcatalog.fits", luminosityBins = n.logspace(38,45,50), Nstack = 400, Nclustering = 400, outputFolder="emissionLineLuminosityFunctions/" , zmin = zmin, zmax = zmax)
		lf.setRedshiftArray( redshiftColumn='Z' )
		lf.setRedshiftSelection( redshiftQualityColumn='ZFLAGS', lowerBound=1.9, upperBound=9.1)
		lf.setWeightArray( 1./(area * lf.catalog['SSR']*lf.catalog['TSR']), area )
		selection = (lf.catalog['TSR']>0) & (lf.catalog['SSR']>0)
		lf.computeHistogramLF(selection)
		lf.computeHistogramVariance(selection,jk=0.1)
		lf.computeMeanWeightedRedshift(selection)
		lf.get_completness_limit(selection)
		lf.correctVolumeEffect(n.ones_like(lf.LF), n.zeros_like(lf.LF))
		lf.writeLF(selection,surveyNameSuffix="WIDE")

import sys
sys.exit()
# from DEEP2 survey
zsDEEP2min=n.array([0.33,0.33,0.4,0.45,0.50,0.60,0.70,0.75,0.78,0.83, 1.16 ])
zsDEEP2max=n.array([0.40,0.45,0.5,0.55,0.60,0.70,0.78,0.8,0.83,1.03, 1.3 ])
linesFittedDEEP2=n.array([ [[O3_5007,"O3_5007"]], [[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"]],[[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"], [O2_mean,"O2_3728"]], [[H1_4862,"H1_4862"],[O2_mean,"O2_3728"]], [[O2_mean,"O2_3728"]], [[O2_mean,"O2_3728"]] ])
#zsDEEP2min=n.array([0.17,0.33,0.33,0.4,0.45,0.50,0.60,0.70,0.78,0.83, 1.16 ])
#zsDEEP2max=n.array([0.36,0.40,0.45,0.5,0.55,0.60,0.70,0.78,0.83,1.03, 1.3 ])
#linesFittedDEEP2=n.array([[[H1_6564,"H1_6564"]], [[O3_5007,"O3_5007"]], [[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"]],[[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[O3_5007,"O3_5007"],[H1_4862,"H1_4862"]], [[H1_4862,"H1_4862"],[O2_mean,"O2_3728"]], [[O2_mean,"O2_3728"]], [[O2_mean,"O2_3728"]] ])
area1=0.60
area2=0.62
area3=0.90
area4=0.66
areaAll=area1+area2+area3+area4

for ii in range(len(zsDEEP2min)):
	zmin = zsDEEP2min[ii]
	zmax = zsDEEP2max[ii]
	lineSet=linesFittedDEEP2[ii]
	for line in lineSet :
		print "---------------------------------------------------"
		print line, zmin, zmax
		lf = LineLuminosityFunction(lineWavelength=line[0], lineName=line[1], cosmology = cosmo, surveyName =  "DEEP2", redshift_catalog = "zcat.deep2.dr4.v2.LFcatalog.fits", luminosityBins = n.logspace(38,45,25), Nstack = 400, Nclustering = 400, outputFolder="emissionLineLuminosityFunctions/" , zmin = zmin, zmax = zmax)
		lf.setRedshiftArray( redshiftColumn='ZBEST' )
		lf.setRedshiftSelection( redshiftQualityColumn='ZQUALITY', lowerBound=0.9, upperBound=7.)
		if zmin < 0.7:
			area = area1
			selection = (lf.catalog['TSR']>0) & (lf.catalog['SSR']>0) & (lf.catalog['DEC']>50.)
			lf.setWeightArray( 1./(area * lf.catalog['SSR']*lf.catalog['TSR']), area )

		if zmin >= 0.7:
			area = areaAll
			lf.setWeightArray( 1./(area * lf.catalog['SSR']*lf.catalog['TSR']), area )
			selection = (lf.catalog['TSR']>0) & (lf.catalog['SSR']>0)

		lf.computeHistogramLF(selection)
		lf.computeHistogramVariance(selection,jk=0.1)
		lf.computeMeanWeightedRedshift(selection)
		lf.get_completness_limit(selection)
		lll=n.array( glob.glob( "/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/" + line[1] + "/" + line[1] + "-GALFORM-DEEP2NEXT-z*Rcorrection.txt" ))

		redshifts=[]
		for el in lll :
			base = el.split('/')[-1]
			out = base.split('-')
			redshifts.append(float(out[3][1:5]))

		redshifts = n.array(redshifts)
		sel = (lf.meanRedshift - 0.05 < redshifts) & (redshifts < lf.meanRedshift + 0.05)
		mCorr, error  = n.loadtxt(lll[sel][0],unpack=True ) 
		lf.correctVolumeEffect(mCorr, error)
		lf.writeLF(selection)




# from VIPERS survey
print "VIPERS"
zsVIMOSmin=n.array([0.65,0.84, 1.1])
zsVIMOSmax=n.array([0.84,1.1, 1.3])
linesFittedVIMOS=n.array([ [[O3_5007,"O3_5007"],[O2_mean,"O2_3728"],[H1_4862,"H1_4862"]], [[O2_mean,"O2_3728"]], [[O2_mean,"O2_3728"]]])

area=24.

for ii in range(len(zsVIMOSmin)):
	zmin = zsVIMOSmin[ii]
	zmax = zsVIMOSmax[ii]
	lineSet=linesFittedVIMOS[ii]
	for line in lineSet :
		lf = LineLuminosityFunction(lineWavelength=line[0], lineName=line[1], cosmology = cosmo, surveyName =  "VIPERS", redshift_catalog = "VIPERS_W14_summary_v1.LFcatalog.fits", luminosityBins = n.logspace(38,45,50), Nstack = 400, Nclustering = 400, outputFolder="emissionLineLuminosityFunctions/" , zmin = zmin, zmax = zmax)
		lf.setRedshiftArray( redshiftColumn='zspec' )
		lf.setRedshiftSelection( redshiftQualityColumn='zflg', lowerBound=0.9, upperBound=100.)
		lf.setWeightArray( 1./(area * lf.catalog['SSR']*lf.catalog['TSR']), area )
		selection = (lf.catalog['TSR']>0) & (lf.catalog['SSR']>0)
		lf.computeHistogramLF(selection)
		print "---------------------------------------------------"
		print line, zmin, zmax, lf.ngals
		lf.computeHistogramVariance(selection,jk=0.1)
		lf.computeMeanWeightedRedshift(selection)
		lf.get_completness_limit(selection)
		lf.correctVolumeEffect(n.ones_like(lf.LF), n.zeros_like(lf.LF))
		lf.writeLF(selection)#,surveyNameSuffix="DEEP")


