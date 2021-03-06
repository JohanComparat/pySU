#! /usr/bin/env python

"""
This script computes the fits the line strength on the spectra from of the MUSE MACS1931 survey
"""

# defines the survey object
from GalaxySurveyMUSE import *
survey = GalaxySurveyMUSE()

# global parameters to add to the survey class
survey.zmin = 0.1 # minimum redshift considered
survey.zmax = 6. # maximum redshift considered
survey.WLmin = 4750. # lambda min
survey.WLmax = 9350. # lambda max of the spectrum
survey.resolution = 3000. # resolution of the Keck deimos setting used
survey.path_to_output_catalog  = join( survey.muse_catalog_dir, "Catalog.spectra_MACS1931.lines.fits")
survey.Ngalaxies=len(survey.catalog)

# line list to be fitted
from lineListAir import *
lineList = n.array([[Ne3,Ne3_3869,"Ne3_3869","left"],[O3,O3_4363,"O3_4363","right"],[O3,O3_4960,"O3_4960","left"],[O3,O3_5007,"O3_5007","right"],[N2,N2_6549,"N2_6549","left"],[N2,N2_6585,"N2_6585","right"],[S2,S2_6718,"S2_6718","left"],[S2,S2_6732,"S2_6732","right"],[Ar3,Ar3_7137,"Ar3_7137","left"]])
recLineList = n.array([[H1,H1_1216,"H1_1216","right"],[H1,H1_3970,"H1_3970","right"],[H1,H1_4102,"H1_4102","right"],[H1,H1_4341,"H1_4341","right"],[H1,H1_4862,"H1_4862","left"],[H1,H1_6564,"H1_6564","left"]])
doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2]])

# import the fitting routines
import LineFittingLibrary as lineFit
lfit  =  lineFit.LineFittingLibrary()

output = n.ones_like(n.empty([survey.Ngalaxies,195]))*lfit.dV
"""
jj=0
catalog_entry = survey.catalog[jj]
spectrum = GalaxySpectrumMUSE(catalog_entry)
spectrum.openObservedSpectrum()
print catalog_entry
print "check flux unit !!"
wl, fl, flErr = spectrum.wavelength, spectrum.fluxl*1e-18, spectrum.fluxlErr*1e-18
d_out,m,h=[],[],[]
datI,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl,flErr,a0= n.array([O2_3727,O2_3729]) *(1+catalog_entry['FINAL_Z']), lineName="O2_3728",p0_sigma=0.1,model="gaussian")
"""
for jj in range(survey.Ngalaxies):
	catalog_entry = survey.catalog[jj]
	spectrum = GalaxySpectrumMUSE(catalog_entry)
	spectrum.openObservedSpectrum()
	if catalog_entry['FINAL_Z']>survey.zmin and catalog_entry['FINAL_Z']<survey.zmax :
		print catalog_entry
		print "check flux unit !!"
		wl, fl, flErr = spectrum.wavelength, spectrum.fluxl*1e-20, spectrum.fluxlErr*1e-20

		if len(wl)>100 : 
			d_out,m,h=[],[],[]
			datI,mI,hI=lfit.fit_Line_OIIdoublet(wl,fl,flErr,a0= n.array([O2_3727,O2_3729]) *(1+catalog_entry['FINAL_Z']), lineName="O2_3728", fitWidth=40, DLC=20,  p0_flux=1e-17, p0_sigma=1.,model="gaussian")
			
			if datI[2] != lfit.dV : 
				lfit.plotLineFit(wl,fl,flErr,mI,n.mean( n.array([O2_3727,O2_3729])* (1+ catalog_entry['FINAL_Z'] ) ),join( survey.muse_catalog_dir,"plots", catalog_entry['SpecName'] + "_O2_3728.png"), title = r"$\log(f)$="+str(n.round(n.log10(datI[2]),2)) + r", $\sigma$="+str(n.round(datI[4],2)) )
			
			d_out.append(datI)
			print len(datI)
			m.append(mI)
			h.append(hI)

			# now fits the emission lines
			for li in lineList:
				datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1+ catalog_entry['FINAL_Z']) , lineName=li[2], continuumSide=li[3], fitWidth=40, DLC=20, p0_flux=2e-16 ,model="gaussian",p0_sigma=1.)
				if datI[1] != lfit.dV : 
					lfit.plotLineFit(wl,fl,flErr,mI,n.mean( n.array([O2_3727,O2_3729])* (1+ catalog_entry['FINAL_Z'] ) ),join( survey.muse_catalog_dir,"plots", catalog_entry['SpecName'] + "_"+li[2]+".png"), title =  r"$\log(f)$="+str(n.round(n.log10(datI[1]),2)) + r", $\sigma$="+str(n.round(datI[3],2)) )
					
				d_out.append(datI)
				m.append(mI)
				h.append(hI)

			# now fits the recombination lines 
			for li in recLineList:
				datI,mI,hI=lfit.fit_Line(wl,fl,flErr,li[1]*(1 +catalog_entry['FINAL_Z'] ),lineName=li[2], continuumSide=li[3], fitWidth=40, DLC=20, p0_flux=2e-16 ,model="gaussian",p0_sigma=1.)
				if datI[1] != lfit.dV : 
					lfit.plotLineFit(wl,fl,flErr,mI,n.mean( n.array([O2_3727,O2_3729])* (1+ catalog_entry['FINAL_Z'] ) ),join( survey.muse_catalog_dir,"plots", catalog_entry['SpecName'] + "_"+li[2]+".png"), title =  r"$\log(abs(f))$="+str(n.round(n.log10(abs(datI[1])),2)) + r", $\sigma$="+str(n.round(datI[3],2)) )

				#print len(datI)
				d_out.append(datI)
				m.append(mI)
				h.append(hI)

			# now outputs the result table
			heading="".join(h)
			out=n.hstack((d_out))
			out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*lfit.dV
			output[jj]=out
			print out

# now writes the output
colNames = heading.split()
hdu2 = fits.BinTableHDU.from_columns(survey.catalog.columns)
new_columns = hdu2.data.columns 

for ii in range(len(colNames)):
	new_columns += fits.Column(name=colNames[ii],format='D', array=output.T[ii] )

hdu = fits.BinTableHDU.from_columns(new_columns)
os.system("rm -rf "+survey.path_to_output_catalog)
hdu.writeto(survey.path_to_output_catalog)
