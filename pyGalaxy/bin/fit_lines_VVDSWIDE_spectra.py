"""
This script computes the fits the line strength on the spectra from of the VVDS DEEP survey
"""
from lib_plot import *
# defines the survey object
import GalaxySurveyVVDS as gs
#survey = gs.GalaxySurveyVVDS(redshift_catalog="VVDS_DEEP_summary.fits")#  
#survey = gs.GalaxySurveyVVDS(redshift_catalog="VVDS_UDEEP_summary.fits")
survey = gs.GalaxySurveyVVDS(redshift_catalog="VVDS_WIDE_summary.fits")

# global parameters to add to the survey class
survey.zmin = 0.01 # minimum redshift considered
survey.zmax = 1.7 # maximum redshift considered
survey.WLmin = 5800. # lambda min
survey.WLmax = 9200. # lambda max of the spectrum
survey.resolution = 200. # resolution of the Keck deimos setting used
#survey.path_to_output_catalog  =  survey.vvds_catalog_dir+ "VVDS_DEEP_summary.linesFitted.fits"
#survey.path_to_output_catalog  =  survey.vvds_catalog_dir+ "VVDS_UDEEP_summary.linesFitted.fits"
survey.path_to_output_catalog  =  survey.vvds_catalog_dir+ "VVDS_WIDE_summary.linesFitted.fits"
survey.Ngalaxies=len(survey.catalog)

plotDir = survey.vvds_dir+'products/emissionLineLuminosityFunctions/plots/lineFits/'

# line list to be fitted
from lineListAir import *
lineList = n.array([ [O2, O2_mean, "O2_3728", "left"], [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"]])

recLineList = n.array([[H1,H1_1216,"H1_1216","right"],[H1,H1_3970,"H1_3970","right"],[H1,H1_4102,"H1_4102","right"],[H1,H1_4341,"H1_4341","right"],[H1,H1_4862,"H1_4862","left"],[H1,H1_6564,"H1_6564","left"]])
#doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])

# import the fitting routines
import LineFittingLibrary as lineFit
lfit  =  lineFit.LineFittingLibrary(fitWidth = 70.)

# import the class to handle spectra
import GalaxySpectrumVVDS as gsp

# construct a magnitude array T05 updated to T07 when possible.
# first loads the filter list
from filterList import *

output = n.ones_like(n.empty([survey.Ngalaxies,194]))*lfit.dV
for jj in range(survey.Ngalaxies):
	catalog_entry = survey.catalog[jj]
	spectrum = gsp.GalaxySpectrumVVDS(catalog_entry,lineFits = True)
	print catalog_entry
	spectrum.openObservedSpectrum()
	wl,fl,flErr = spectrum.wavelength, spectrum.fluxl, spectrum.fluxlErr

	if catalog_entry['Z']>survey.zmin and catalog_entry['Z']<survey.zmax and len(wl)>200:
		# corrects for the slit aperture to the magnitude observed
		foI,foI_err = lfit.getFractionObsMag(mag=catalog_entry['MAGI'], lambdaMag=lambIcfht, filter=filterIcfht, xmin=n.max([5800.,n.min(wl)]), xmax=n.min([9200.,n.max(wl)]), wl=wl,fl=fl,flErr=flErr)

		#initialises the array where the fits are stored
		d_out,m,h=[],[],[]
		"""
		# now fits the emission lines
		datI,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl,flErr,a0= n.array([O2_3727,O2_3729])* (1+catalog_entry['Z']),lineName="O2_3728",p0_sigma=1,model="gaussian")
		d_out.append(datI)
		m.append(mI)
		h.append(hI)
		"""
		for li in lineList:
			datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1+catalog_entry['Z']), lineName=li[2], continuumSide=li[3], model="gaussian",p0_sigma=1)
			d_out.append(datI)
			m.append(mI)
			h.append(hI)
			"""
			if datI[1]>0 and (li[2]=="O2_3728" or li[2]=="O3_5007"):
				p.figure(0,(5,5))
				p.axes([0.2,0.23,0.75,0.72])
				p.errorbar(wl,fl,yerr=flErr,label='data')
				p.plot(mI[0],mI[1],label='model')
				p.xlim((mI[0].min()-10,mI[0].max()+10))
				p.ylim((mI[1].min()/5,mI[1].max()*5))
				p.xticks(rotation=60)
				p.xlabel('Wavelength [A]')
				p.ylabel('Flux density')
				p.yscale('log')
				p.legend(fontsize=9)
				p.title(str(datI[1]) +", "+ str(datI[2]),fontsize=9)
				p.savefig(plotDir+'spec_'+str(catalog_entry['NUM'])+'_'+li[2]+'.pdf')
				p.clf()
			"""

		# now fits the recombination lines 
		for li in recLineList:
			datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1+catalog_entry['Z']), lineName=li[2], continuumSide=li[3],model="gaussian",p0_sigma=1)
			d_out.append(datI)
			m.append(mI)
			h.append(hI)

		# now outputs the result table
		head=n.hstack((["fo fo_err"],h))
		heading="".join(head)
		out=n.hstack((n.hstack((n.array([foI, foI_err]))), n.hstack((d_out))))
		out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*lfit.dV
		output[jj]=out
		#print out

# now writes the output
colNames = heading.split()
hdu2 = gs.fits.BinTableHDU.from_columns(survey.catalog.columns)
new_columns = hdu2.data.columns 

for jj in range(len(colNames)):
	new_columns += gs.fits.Column(name=colNames[jj],format='D', array=output.T[jj] )

hdu = gs.fits.BinTableHDU.from_columns(new_columns)
hdu.writeto(survey.path_to_output_catalog)
