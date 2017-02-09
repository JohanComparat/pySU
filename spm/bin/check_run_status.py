
#plate = plates[0]
tbdata = hdus[1].data

for plate in plates:
	print plate, time.time()
	mask = tbdata['PLATE'] ==int(plate)
	if len(mask.nonzero()[0])>=1 :
		newtbdata = tbdata[mask]
		hdu = fits.BinTableHDU(data=newtbdata)
		newTab=join( os.environ['SDSSDR12_DIR'], "catalogs", "perPlate", "sp-"+plate+".fits")
		hdu.writeto(newTab)

from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits

plate = sys.argv[1]
dir ='stellarpop-m11-salpeter'

hdus = fits.open( join( os.environ['SDSSDR12_DIR'], "catalogs", "specObj-dr12.fits") )

#plates_all = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "plateNumberList"), unpack=True, dtype='str')
#plates = plates_all[:-2]


dV=-9999.99
def get_table_entry_full(hduSPM):
	headerA =" age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
	
	table_entry = [ 10**hduSPM.header['age_universe'], 10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean_up']-10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean']-10**hduSPM.header['age_lightW_mean_low'], hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean_up'] - hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean'] - hduSPM.header['metallicity_lightW_mean_low'], 10**hduSPM.header['age_massW_mean'], 10**hduSPM.header['age_massW_mean_up']-10**hduSPM.header['age_massW_mean'], 10**hduSPM.header['age_massW_mean']-10**hduSPM.header['age_massW_mean_low'], hduSPM.header['metallicity_massW_mean'], hduSPM.header['metallicity_massW_mean_up'] - hduSPM.header['metallicity_massW_mean'], hduSPM.header['metallicity_massW_mean'] - hduSPM.header['metallicity_massW_mean_low'], hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean_up'] - hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean'] - hduSPM.header['stellar_mass_mean_low'], hduSPM.header['EBV'], hduSPM.header['ssp_number']]
	#print hduSPM.header
	for iii in n.arange(hduSPM.header['ssp_number']):
		table_entry.append( hduSPM.header['stellar_mass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['age_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['metal_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['SFR_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightMass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightLight_ssp_'+str(iii)] )
		headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)
	
	if hduSPM.header['ssp_number']<8 :
		for iii in n.arange(hduSPM.header['ssp_number'], 8, 1):
			table_entry.append([dV, dV, dV, dV, dV, dV])
			headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

	table_entry = n.array( n.hstack((table_entry)) )
	#print table_entry.shape
	return n.hstack((table_entry)), headerA
	
# step 2 : match to thecreated data set
init_cat = join( os.environ['SDSSDR12_DIR'], "catalogs", "perPlate", "sp-"+plate+".fits")
plate_catalog = join( os.environ['SDSSDR12_DIR'], dir, "catalogs", "spFly-"+plate+".fits")
	
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

table_all = []

for fiber, mjd in zip(orig_table['FIBERID'], orig_table['MJD']):
	fitFile = join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", plate, "spec-"+plate+"-"+str(mjd)+"-"+str(fiber).zfill(4)+"-SPM-MILES.fits")
	if os.path.isfile(fitFile):
		table_entry, headers = get_table_entry_full( hduSPM=fits.open(fitFile)[1] )
		table_all.append(table_entry)
	else:
		table_all.append(n.ones(66)*dV)

newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto(plate_catalog)


def make_summary_table(plate, dir = dir):

plate=plates[-1]
t0 = time.time()
data = hdus[1].data[(hdus[1].data['PLATE']==int(plate))]
print time.time()-t0 #30 seconde
t1 = time.time()
fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", plate, '*.fits')))

fitFile = fitList[0]
mjd, fiber = os.path.basename(fitFile).split('-')[2:4]
hd = fits.open(fitFile)
newcols=[]

for head in hd[1].header.keys()[14:]:
	print head
	data.columns.add_col( fits.Column(name=head, format="D", array = n.ones_like(data['PLATE']) * -99.99) )


cols = fits.ColDefs(data.columns)
tbhdu = fits.BinTableHDU.from_columns(cols)

sel=(data['MJD']==int(mjd)) & (data['FIBERID']==int(fiber))

data[head][sel] = hd[1].header[head]

def get_model_result(fitFile):



modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", plate, '*.model')))
print time.time()-t1 #30 seconde




def get_lists_fits_models(plate, dir = dir):
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(int(plate)), '*.fits')))
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(int(plate)), '*.model')))
	return len(fitList), len(modList)

def get_lists_fits_tables(plate, dir = dir):
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "tables", str(int(plate))+ '_full.data')))
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(int(plate)), '*.model')))
	return len(fitList), len(modList)

"""
nF=n.zeros_like(plates)
nM=n.zeros_like(plates)
for ii, plate in enumerate(plates):
	nF[ii], nM[ii] = get_lists_fits_models(plate, dir=dir)
	#print ii, plate, nF[ii], nM[ii]#, nM[ii] *100./nF[ii]

plates[(nM<nF)]

 
nF=n.zeros_like(plates)
nM=n.zeros_like(plates)
for ii, plate in enumerate(plates):
	nF[ii], nM[ii] = get_lists_fits_tables(plate, dir=dir)
	print ii, plate, nF[ii], nM[ii]#, nM[ii] *100./nF[ii]

plates[(nF=='0')]

n.savetxt("/users/comparat/batchscripts_firefly_salpeter_table/plates_to_run.ascii", n.transpose(plates[(nF=='0')]), fmt='%s')
"""
def get_plate_lists(plate, dir = dir):
	specList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], 'spectra', str(plate), 'spec-*.fits')))
	specList.sort()
	
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(plate), '*.fits')))
	fitList.sort()	
	
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(plate), '*.model')))
	modList.sort()	
	
	tabList_l = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "tables", str(plate)+'*line*.data')))
	tabList_l.sort()	
	
	tabList_f = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "tables", str(plate)+'*full*.data')))
	tabList_f.sort()
	
	return specList, fitList, modList, tabList_l, tabList_f


def get_unprocessed_fiber_lists_per_plate(plate, dir = dir):
	in_plate = (hdus[1].data['PLATE']==int(plate))
	gal = (in_plate) & (hdus[1].data['CLASS']=="GALAXY") & (hdus[1].data['Z']>0) & (hdus[1].data['Z']<1.7)
	fiber_2_fit = hdus[1].data['FIBERID'][in_plate]
	mjd_2_fit = hdus[1].data['MJD'][in_plate]
	
	specList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], 'spectra', str(plate), 'spec-*.fits')))
	specList.sort()
	fib_spec = n.array([int(os.path.basename(fl).split('-')[3][:-5]) for fl in specList ])
	mjd_spec = n.array([int(os.path.basename(fl).split('-')[2]) for fl in specList ])
	
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(plate)+'*.fits')))
	fitList.sort()	
	fib_fitted = n.array([int(os.path.basename(fl).split('-')[3][:-5]) for fl in fitList ])
	mjd_fitted = n.array([int(os.path.basename(fl).split('-')[2]) for fl in fitList ])
	remaining_fibers = set(fib_fitted).difference(set(fiber_2_fit))
	return remaining_fibers 

def get_info_from_catalog(plate):
	in_plate = (hdus[1].data['PLATE']==int(plate))
	gal = (in_plate) & (hdus[1].data['ZWARNING']==0) & (hdus[1].data['CLASS_NOQSO']=="GALAXY") & (hdus[1].data['Z_NOQSO'] > hdus[1].data['Z_ERR_NOQSO']) & (hdus[1].data['Z_ERR_NOQSO']>0)
	#spec.hdulist[2].data['CLASS_NOQSO'][0]=="GALAXY" and spec.hdulist[2].data['Z_NOQSO'][0] >  spec.hdulist[2].data['Z_ERR_NOQSO'][0] and spec.hdulist[2].data['Z_ERR_NOQSO'][0]>0 and spec.hdulist[2].data['ZWARNING'][0] ==0 and abs(spec.hdulist[2].data['Z_NOQSO'][0] - spec.hdulist[2].data['Z'][0])>abs(spec.hdulist[2].data['Z_ERR_NOQSO'][0])
	return len(hdus[1].data['Z'][gal])

def get_plate_lists_light(plate, dir = dir):
	specList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], 'spectra', str(plate), 'spec-*.fits')))
	specList.sort()
	
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(plate), '*.fits')))
	fitList.sort()	
	
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(plate), '*.model')))
	modList.sort()	
	
	return specList, fitList, modList

def create_light_table(plates, outname = "run-status-table.ascii",  dir = dir):
	Nspec = n.zeros(len(plates))
	Nmodel = n.zeros(len(plates))
	Nfit = n.zeros(len(plates))
	for ii, plate in enumerate(plates):
		specList, fitList, modList = get_plate_lists_light(plate)
		Nspec[ii]=len(specList)
		Nfit[ii]=len(fitList)
		Nmodel[ii]=len(modList)
		
	n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, outname), n.transpose([plates.astype('int'), Nspec, Nfit,Nmodel]), header='plate Nspec Nfit Nmodel')
			
def create_basic_table(plates, outname = "run-status-table.ascii",  dir = dir):
	Nspec = n.zeros(len(plates))
	Ngal = n.zeros(len(plates))
	Nmodel = n.zeros(len(plates))
	Nfit = n.zeros(len(plates))
	lenTableLine = n.zeros(len(plates))
	lenTableFull = n.zeros(len(plates))
	for ii, plate in enumerate(plates):
		t0 = time.time()
		specList, fitList, modList, tabList_l, tabList_f = get_plate_lists(plate)
		n_fittable_gal = get_info_from_catalog(plate)
		Nspec[ii]=len(specList)
		Nfit[ii]=len(fitList)
		Ngal[ii]= n_fittable_gal
		Nmodel[ii]=len(modList)
		if len(tabList_l)>0:
			tab_l = n.loadtxt(tabList_l[0],unpack=True)
			tab_f = n.loadtxt(tabList_f[0],unpack=True)
			lenTableLine[ii] = len(tab_l)
			lenTableFull[ii] = len(tab_f)
		print ii, plate, time.time() - t0
	n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, outname), n.transpose([plates.astype('int'), Nspec, Ngal, Nfit, Nmodel, lenTableLine, lenTableFull]), header='plate Nspec Ngal Nfit Nmodel NinTable_l NinTable_f')
			
create_basic_table(plates)

# create_light_table(plates, "run-status-table-light.ascii")
	
for plate in plates:
	print plate
	fibs = get_unprocessed_fiber_lists_per_plate(plate)
	if len(fibs)>0:
		n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, "remaining_2_process", "remaining_fibers_"+plate+".txt"), n.transpose([n.ones_like(n.array(list(fibs)))*int(plate),n.array(list(fibs))]))
	
# now exploit the data created
outname = "run-status-table.ascii"
plates, Nspec, Ngal, Nfit, Nmodel, lenTableLine,lenTableFull = n.loadtxt(join( os.environ['SDSSDR12_DIR'], dir, outname), unpack= True)

comp = 0.9
isSpec = (Nspec > 0)
isGal = (Ngal > 0) & (isSpec)
not_fitted_FF = ( Nfit < Ngal ) & (isGal)
done_model = (Ngal * comp <= Nmodel) & (isGal)
done_table_l = (Ngal * comp <= lenTableLine) & (isGal)
done_table_f = (Ngal * comp <= lenTableFull) & (isGal)

f=open(join( os.environ['SDSSDR12_DIR'], dir, "status-30-12.txt"), 'w')
f.write("Status of the run \n")
f.write("----------------------------------------\n")
f.write("Author: JC \n")
f.write("last update: "+time.ctime()+" \n")
f.write("----------------------------------------\n")
f.write("Total data: " + str(len(Nspec[isSpec])) + " plates containing " + str(int(n.sum(Nspec)))+ " spectra (contains sky and std star fibers) \n")
f.write("Total galaxy: "+ str(int(n.sum(Ngal))) + " galaxies with CLASS=GALAXY AND 0<z<1.7 \n")
f.write("Total FFfit: " + str(int(n.sum(Nfit))) + " , i.e. "+str(n.round(100.*n.sum(Nfit)/n.sum(Ngal),2))+" per cent \n")

# incomplete plate list
print len(plates[not_fitted_FF])

f.write("Total emission line: " + str(int(n.sum(Nmodel)))+" spectra have an emission line model, i.e. "+str(n.round(100.*n.sum(Nmodel)/n.sum(Ngal),2))+" per cent \n")

f.write( str(int(n.sum(lenTableLine)))+" spectra are in a summary table, i.e. "+str(n.round(100.*n.sum(lenTableLine)/n.sum(Ngal),2))+" per cent \n")

f.close()



sys.exit()

def runSpec(specLiteFile, dir = dir):
	baseN = os.path.basename(specLiteFile).split('-')
	plate = baseN[1] #7457# sys.argv[1] #7619
	mjd = baseN[2] # 56746#sys.argv[2] # 56900
	fibre = baseN[3] # 471#sys.argv[3] #300
	outputFolder = join( os.environ['SDSSDR12_DIR'], dir, 'stellarpop', str(plate))
	if os.path.isdir(outputFolder)==False:
		os.mkdir(outputFolder)

for plate in plates:
	rootName = join(os.environ['HOME'], "batchscripts_firefly", plate)
	writeScript(rootName, plate)
	

for el in fileList:
	spec = runSpec(el)