from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits
from os.path import join
import astropy.cosmology as co
cosmo = co.Planck13
import astropy.io.fits as fits

# for one galaxy spectrum
import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm

init_cat=join(os.environ['DEEP2_DIR'], "catalogs", "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.fits")
summary_catalog = join(os.environ['DEEP2_DIR'], "catalogs", "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.spm.fits")
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

kroupaFolder = join( os.environ['DEEP2_DIR'], 'stellarpop-m11-kroupa', 'stellarpop')
salpeterFolder = join( os.environ['DEEP2_DIR'], 'stellarpop-m11-salpeter', 'stellarpop')


dV=-9999.99

def get_table_entry_full(hduSPM):
	# print "gets entry"
	hduSPM.header
	prefix = hduSPM.header['IMF'] + "_" + hduSPM.header['library'] + "_"
	#print prefix
	headerA =" "+prefix+"age_lightW "+prefix+"age_lightW_err_plus "+prefix+"age_lightW_err_minus "+prefix+"metallicity_lightW "+prefix+"metallicity_lightW_err_plus "+prefix+"metallicity_lightW_err_minus "+prefix+"age_massW "+prefix+"age_massW_err_plus "+prefix+"age_massW_err_minus "+prefix+"metallicity_massW "+prefix+"metallicity_massW_err_plus "+prefix+"metallicity_massW_err_minus "+prefix+"stellar_mass "+prefix+"stellar_mass_err_plus "+prefix+"stellar_mass_err_minus "+prefix+"spm_EBV "+prefix+"nComponentsSSP "+prefix+"chi2 "+prefix+"ndof "
	
	table_entry = [10**hduSPM.header['age_lightW']          
	, 10**hduSPM.header['age_lightW_up']       
	, 10**hduSPM.header['age_lightW_low']      
	, hduSPM.header['metallicity_lightW']  
	, hduSPM.header['metallicity_lightW_up']
	, hduSPM.header['metallicity_lightW_low']
	, 10**hduSPM.header['age_massW']           
	, 10**hduSPM.header['age_massW_up']        
	, 10**hduSPM.header['age_massW_low']       
	, hduSPM.header['metallicity_massW']   
	, hduSPM.header['metallicity_massW_up']
	, hduSPM.header['metallicity_massW_low']
	, hduSPM.header['stellar_mass']        
	, hduSPM.header['stellar_mass_up']     
	, hduSPM.header['stellar_mass_low']    
	, hduSPM.header['EBV'] 
	, hduSPM.header['ssp_number']
	, hduSPM.header['chi2']
	, hduSPM.header['ndof']
	]
	
	#print hduSPM.header
	for iii in n.arange(hduSPM.header['ssp_number']):
		table_entry.append( hduSPM.header['stellar_mass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['age_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['metal_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['SFR_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightMass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightLight_ssp_'+str(iii)] )
		headerA += ' '+prefix+'stellar_mass_ssp_'+str(iii) + ' '+prefix+'ndofage_ssp_'+str(iii) + ' '+prefix+'metal_ssp_'+str(iii) + ' '+prefix+'SFR_ssp_'+str(iii) + ' '+prefix+'weightMass_ssp_'+str(iii) + ' '+prefix+'weightLight_ssp_'+str(iii)
	
	if hduSPM.header['ssp_number']<8 :
		for iii in n.arange(hduSPM.header['ssp_number'], 8, 1):
			table_entry.append([dV, dV, dV, dV, dV, dV])
			headerA += ' '+prefix+'stellar_mass_ssp_'+str(iii) + ' '+prefix+'age_ssp_'+str(iii) + ' '+prefix+'metal_ssp_'+str(iii) + ' '+prefix+'SFR_ssp_'+str(iii) + ' '+prefix+'weightMass_ssp_'+str(iii) + ' '+prefix+'weightLight_ssp_'+str(iii)

	table_entry = n.array( n.hstack((table_entry)) )
	#print table_entry.shape
	return n.hstack((table_entry)), headerA
	


table_all = []
table_entry = []
headers = ""
headers = ""
for mask, objno in zip(orig_table['MASK'], orig_table['OBJNO']):
	fitFile = join( os.environ['DEEP2_DIR'], 'stellarpop', "spFly-deep2-"+str(mask)+"-"+str(objno)+".fits")
	if os.path.isfile(fitFile):
		print fitFile
		hdus = fits.open(fitFile)
		for ii in range(1,len(hdus)):
			table_entry_i, headers_i = get_table_entry_full( hdus[ii] )
			table_entry.append(table_entry_i)
			headers += headers_i
		
		#print len(n.hstack((table_entry)))
		table_all.append(n.hstack((table_entry)))
		#print headers
		#print table_all[-1]
		#print len(table_all[-1])
		#fitFileLast = fitFile
	else:
		table_all.append(n.ones(603)*dV)

newDat = n.transpose(table_all)
all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
if os.path.isfile(summary_catalog):
	os.remove(summary_catalog)

hdu.writeto(summary_catalog)

