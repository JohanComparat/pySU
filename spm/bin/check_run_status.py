from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits
dir ='stellarpop-m11-salpeter'
hdus = fits.open( join( os.environ['SDSSDR12_DIR'], "catalogs", "specObj-dr12.fits") )

def get_lists_fits_models(plate, dir ='stellarpop-m11-salpeter'):
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(int(plate)), '*.fits')))
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(int(plate)), '*.model')))
	return len(fitList), len(modList)

"""
nF=n.zeros_like(plates)
nM=n.zeros_like(plates)
for ii, plate in enumerate(plates):
	nF[ii], nM[ii] = get_lists_fits_models(plate, dir=dir)
	print ii, plate, nF[ii], nM[ii], nM[ii] *100./nF[ii]

"""
def get_plate_lists(plate, dir ='stellarpop-m11-salpeter'):
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


def get_unprocessed_fiber_lists_per_plate(plate, dir ='stellarpop-m11-salpeter'):
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
	gal = (in_plate) & (hdus[1].data['CLASS']=="GALAXY") & (hdus[1].data['Z']>0) & (hdus[1].data['Z']<1.7)
	return len(hdus[1].data['Z'][gal])

def get_plate_lists_light(plate, dir ='stellarpop-m11-salpeter'):
	specList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], 'spectra', str(plate), 'spec-*.fits')))
	specList.sort()
	
	fitList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", str(plate), '*.fits')))
	fitList.sort()	
	
	modList = n.array(glob.glob(join( os.environ['SDSSDR12_DIR'], dir, "model", str(plate), '*.model')))
	modList.sort()	
	
	return specList, fitList, modList

def create_light_table(plates, outname = "run-status-table.ascii",  dir ='stellarpop-m11-salpeter'):
	Nspec = n.zeros(len(plates))
	Nmodel = n.zeros(len(plates))
	Nfit = n.zeros(len(plates))
	for ii, plate in enumerate(plates):
		specList, fitList, modList = get_plate_lists_light(plate)
		Nspec[ii]=len(specList)
		Nfit[ii]=len(fitList)
		Nmodel[ii]=len(modList)
		
	n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, outname), n.transpose([plates.astype('int'), Nspec, Nfit,Nmodel]), header='plate Nspec Nfit Nmodel')
			
def create_basic_table(plates, outname = "run-status-table.ascii",  dir ='stellarpop-m11-salpeter'):
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
			
plates_all = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "plateNumberList"), unpack=True, dtype='str')
plates = plates_all[:-2]

#create_basic_table(plates)

# create_light_table(plates, "run-status-table-light.ascii")
	
for plate in plates:
	print plate
	fibs = get_unprocessed_fiber_lists_per_plate(plate)
	if len(fibs)>0:
		n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, "remaining_2_process", "remaining_fibers_"+plate+".txt"), n.transpose([n.ones_like(n.array(list(fibs)))*int(plate),n.array(list(fibs))]))
	
# now exploit the data created
outname = "run-status-table.ascii"
dir ='stellarpop-m11-salpeter'
plates, Nspec, Ngal, Nfit, Nmodel, lenTableLine,lenTableFull = n.loadtxt(join( os.environ['SDSSDR12_DIR'], dir, outname), unpack= True)

comp = 0.9
isSpec = (Nspec > 0)
isGal = (Ngal > 0) & (isSpec)
not_fitted_FF = ( Nfit < Ngal ) & (isGal)
done_model = (Ngal * comp <= Nmodel) & (isGal)
done_table_l = (Ngal * comp <= lenTableLine) & (isGal)
done_table_f = (Ngal * comp <= lenTableFull) & (isGal)

f=open(join( os.environ['SDSSDR12_DIR'], dir, "status.txt"), 'w')
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

def runSpec(specLiteFile):
	baseN = os.path.basename(specLiteFile).split('-')
	plate = baseN[1] #7457# sys.argv[1] #7619
	mjd = baseN[2] # 56746#sys.argv[2] # 56900
	fibre = baseN[3] # 471#sys.argv[3] #300
	outputFolder = join( os.environ['SDSSDR12_DIR'], 'stellarpop-m11-salpeter', 'stellarpop', str(plate))
	if os.path.isdir(outputFolder)==False:
		os.mkdir(outputFolder)

for plate in plates:
	rootName = join(os.environ['HOME'], "batchscripts_firefly", plate)
	writeScript(rootName, plate)
	

for el in fileList:
	spec = runSpec(el)