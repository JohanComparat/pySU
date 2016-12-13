from os.path import join
import os
import numpy as n
import glob 
import sys 

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
	Nmodel = n.zeros(len(plates))
	Nfit = n.zeros(len(plates))
	lenTableLine = n.zeros(len(plates))
	lenTableFull = n.zeros(len(plates))
	for ii, plate in enumerate(plates):
		specList, fitList, modList, tabList_l, tabList_f = get_plate_lists(plate)
		Nspec[ii]=len(specList)
		Nfit[ii]=len(fitList)
		Nmodel[ii]=len(modList)
		if len(tabList_l)>0:
			tab_l = n.loadtxt(tabList_l[0],unpack=True)
			tab_f = n.loadtxt(tabList_f[0],unpack=True)
			lenTableLine[ii] = len(tab_l)
			lenTableFull[ii] = len(tab_f)
	n.savetxt(join( os.environ['SDSSDR12_DIR'], dir, outname), n.transpose([plates.astype('int'), Nspec, Nfit,Nmodel, lenTableLine, lenTableFull]), header='plate Nspec Nfit Nmodel NinTable_l NinTable_f')
			
plates_all = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "plateNumberList"), unpack=True, dtype='str')
plates = plates_all[:-2]
#create_basic_table(plates)
create_light_table("run-status-table-light.ascii")

# now exploit the data created
outname = "run-status-table.ascii"
dir ='stellarpop-m11-salpeter'
plates, Nspec, Nfit, Nmodel, lenTableLine,lenTableFull = n.loadtxt(join( os.environ['SDSSDR12_DIR'], dir, outname))

isSpec = (Nspec > 0)
done_fit = (Nspec == Nfit) & (isSpec)
done_model = (Nspec == Nmodel) & (isSpec)
done_table_l = (Nspec == lenTableLine) & (isSpec)
done_table_f = (Nspec == lenTableFull) & (isSpec)

print len(Nspec[isSpec]), "plates for fit conatining ", n.sum(Nspec)," spectra."
print len(Nspec[done_fit]), "plates were fitted by firefly."
print len(Nspec[done_model]), "plates were line-modeled."
print len(Nspec[done_table_l]), "plates have a line table."
print len(Nspec[done_table_f]), "plates have a full table."

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