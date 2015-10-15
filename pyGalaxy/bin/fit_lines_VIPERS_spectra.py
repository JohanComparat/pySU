"""
This script computes the fits the line strength on the spectra from of the VIPERS survey
"""

# defines the survey object
import GalaxySurveyVIPERS as gs
survey = gs.GalaxySurveyVIPERS(redshift_catalog="VIPERS_W14_summary_v1.fits")

# global parameters to add to the survey class
survey.zmin = 0.5 # minimum redshift considered
survey.zmax = 1.4 # maximum redshift considered
survey.WLmin = 5800. # lambda min
survey.WLmax = 9200. # lambda max of the spectrum
survey.resolution = 200. # resolution of the Keck deimos setting used
survey.path_to_output_catalog  =  survey.vipers_catalog_dir+ "VIPERS_W14_summary_v1.linesFitted.fits"
survey.Ngalaxies=len(survey.catalog)

# line list to be fitted
from lineListAir import *
lineList = n.array([ [O2, O2_mean, "O2_3728", "left"], [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"]])
recLineList = n.array([[H1,H1_1216,"H1_1216","right"],[H1,H1_3970,"H1_3970","right"],[H1,H1_4102,"H1_4102","right"],[H1,H1_4341,"H1_4341","right"],[H1,H1_4862,"H1_4862","left"],[H1,H1_6564,"H1_6564","left"]])
#doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])

# import the fitting routines
import LineFittingLibrary as lineFit
lfit  =  lineFit.LineFittingLibrary(fitWidth = 70.)

# import the class to handle spectra
import GalaxySpectrumVIPERS as gsp

# construct a magnitude array T05 updated to T07 when possible.
# first loads the filter list
from filterList import *

mag0=n.zeros(survey.Ngalaxies)
mag1=n.zeros(survey.Ngalaxies)
mag2=n.zeros(survey.Ngalaxies)
mag3=n.zeros(survey.Ngalaxies)
mag4=n.zeros(survey.Ngalaxies)

mCFHT07=(survey.catalog['u_T07']>12.)&(survey.catalog['u_T07']<29.)
mCFHT05=(survey.catalog['u']>12.)&(survey.catalog['u']<29.)&(mCFHT07==False)
mNoPhot=(mCFHT05==False)&(mCFHT07==False)
mag0[mCFHT07]=survey.catalog['u_T07'][mCFHT07]
mag0[mCFHT05]=survey.catalog['u'][mCFHT05]
mag0[mNoPhot]=n.ones_like(survey.catalog['u'][mNoPhot])*lfit.dV

mCFHT07=(survey.catalog['g_T07']>12.)&(survey.catalog['g_T07']<29.)
mCFHT05=(survey.catalog['g']>12.)&(survey.catalog['g']<29.)&(mCFHT07==False)
mNoPhot=(mCFHT05==False)&(mCFHT07==False)
mag1[mCFHT07]=survey.catalog['g_T07'][mCFHT07]
mag1[mCFHT05]=survey.catalog['g'][mCFHT05]
mag1[mNoPhot]=n.ones_like(survey.catalog['g'][mNoPhot])*lfit.dV

mCFHT07=(survey.catalog['r_T07']>12.)&(survey.catalog['r_T07']<29.)
mCFHT05=(survey.catalog['r']>12.)&(survey.catalog['r']<29.)&(mCFHT07==False)
mNoPhot=(mCFHT05==False)&(mCFHT07==False)
mag2[mCFHT07]=survey.catalog['r_T07'][mCFHT07]
mag2[mCFHT05]=survey.catalog['r'][mCFHT05]
mag2[mNoPhot]=n.ones_like(survey.catalog['r'][mNoPhot])*lfit.dV

mCFHT07=(survey.catalog['i_T07']>12.)&(survey.catalog['i_T07']<29.)
mCFHT05=(survey.catalog['i']>12.)&(survey.catalog['i']<29.)&(mCFHT07==False)
mag3[mCFHT07]=survey.catalog['i_T07'][mCFHT07]
mag3[mCFHT05]=survey.catalog['i'][mCFHT05]

mCFHT07=(survey.catalog['z_T07']>12.)&(survey.catalog['z_T07']<29.)
mCFHT05=(survey.catalog['z']>12.)&(survey.catalog['z']<29.)&(mCFHT07==False)
mNoPhot=(mCFHT05==False)&(mCFHT07==False)
mag4[mCFHT07]=survey.catalog['z_T07'][mCFHT07]
mag4[mCFHT05]=survey.catalog['z'][mCFHT05]
mag4[mNoPhot]=n.ones_like(survey.catalog['z'][mNoPhot])*lfit.dV


#sel=(survey.catalog['NUM']==112154559)
#catalog_entry = survey.catalog[sel][0]
output = n.ones_like(n.empty([survey.Ngalaxies,203]))*lfit.dV
for jj in range(survey.Ngalaxies):
	catalog_entry = survey.catalog[jj]
	spectrum = gsp.GalaxySpectrumVIPERS(catalog_entry,calibration = False,lineFits = True)
	print catalog_entry
	spectrum.openObservedSpectrum()
	wl,fl,flErr = spectrum.wavelength, spectrum.fluxl, spectrum.fluxlErr

	if catalog_entry['zspec']>survey.zmin and catalog_entry['zspec']<survey.zmax and len(wl)>200:
		# corrects for the slit aperture to the magnitude observed
		foR,foR_err = lfit.getFractionObsMag(mag=mag2[jj], lambdaMag=lambRcfht, filter=filterRcfht, xmin=n.max([5800.,n.min(wl)]), xmax=n.min([9200.,n.max(wl)]), wl=wl,fl=fl,flErr=flErr)
		foI,foI_err = lfit.getFractionObsMag(mag=mag3[jj], lambdaMag=lambIcfht, filter=filterIcfht, xmin=n.max([5800.,n.min(wl)]), xmax=n.min([9200.,n.max(wl)]), wl=wl,fl=fl,flErr=flErr)
		foIndex=n.argmax([ foR/foR_err, foI/foI_err])
		fos=n.array([ foR, foI])
		foErrs=n.array([ foR_err, foI_err])
		fo=fos[foIndex]
		fo_err=foErrs[foIndex]

		#initialises the array where the fits are stored
		d_out,m,h=[],[],[]
		# now fits the emission lines
		#datI,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl,flErr,a0= n.array([O2_3727,O2_3729])* (1+catalog_entry['zspec']),lineName="O2_3728",p0_sigma=1,model="gaussian")
		#d_out.append(datI)
		#m.append(mI)
		#h.append(hI)
		for li in lineList:
			datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1+catalog_entry['zspec']), lineName=li[2], continuumSide=li[3],model="gaussian",p0_sigma=1)
			d_out.append(datI)
			m.append(mI)
			h.append(hI)
		# now fits the recombination lines 
		for li in recLineList:
			datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1+catalog_entry['zspec']), lineName=li[2], continuumSide=li[3],model="gaussian",p0_sigma=1)
			d_out.append(datI)
			m.append(mI)
			h.append(hI)

		# now outputs the result table
		head=n.hstack((["foR foR_err foI foI_err fo fo_err MAGU MAGG MAGR MAGI MAGZ"],h))
		heading="".join(head)
		out=n.hstack((n.hstack((n.array([foR, foR_err, foI, foI_err, fo, fo_err, mag0[jj], mag1[jj], mag2[jj], mag3[jj], mag4[jj]]))), n.hstack((d_out))))
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
