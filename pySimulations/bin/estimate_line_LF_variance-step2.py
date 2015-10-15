"""
This script produces quality plots to check that the LFs are fine.
"""

import sys
from lib_plot import *
from lineListAir import *
import glob


plotDir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/plots/"
dir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/*/"
lf_measurement_files_means = n.array(glob.glob(dir+"*NEXT-z*.txt"))
lf_measurement_files_realizations = n.array(glob.glob(dir+"*NEXT-R*.txt"))

elements,surveys,redshifts=[],[],[]
for el in lf_measurement_files_means :
	base = el.split('/')[-1]
	out = base.split('-')
	elements.append(out[0])
	surveys.append(out[2])
	redshifts.append(float(out[3][1:-4]))

elements = n.array(elements)
surveys = n.array(surveys)
redshifts = n.array(redshifts)


element,survey,RR,redshift=[],[],[],[]
for el in lf_measurement_files_realizations :
	base = el.split('/')[-1]
	out = base.split('-')
	element.append(out[0])
	survey.append(out[2])
	RR.append(int(out[3][1:]))
	redshift.append(float(out[4][1:-4]))

element = n.array(element)
survey = n.array(survey)
RR = n.array(RR)
redshift = n.array(redshift)

for jj in range(len(elements)):
	print lf_measurement_files_means[jj]
	selection=(survey == surveys[jj]) & (element == elements[jj]) & (redshift> redshifts[jj] - 0.02) & (redshift < redshifts[jj] + 0.02)
	phis = []
	for el in lf_measurement_files_realizations[selection] :
		Lmin, Lmax, Lmean, phi, phi_err_poisson, phi_err_jackknife, ngals = n.loadtxt(el,unpack=True)
		phis.append(phi)

	phis=n.array(phis)
	phiVar = n.std(phis,axis=0) # variance on phi when using a small area
	phisMean = n.mean(phis,axis=0) # mean phi obtained if small area
	Lmin, Lmax, Lmean, phi, phi_err_poisson, phi_err_jackknife, ngals = n.loadtxt(lf_measurement_files_means[jj], unpack=True) # value using the complete volume
	n.savetxt(lf_measurement_files_means[jj][:-4]+"Rcorrection.txt", n.transpose([ phi / phisMean, phiVar / phisMean ]), header = " multiplicativeCorrection errorPercentageFromSimulation " )

