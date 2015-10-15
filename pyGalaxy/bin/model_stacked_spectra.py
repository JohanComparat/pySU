from ModelSpectraStacks import *
import glob
stack_files=n.hstack(( n.array(glob.glob( "../../LFstacks/*DEEP2*linesFittedfireflyFitsMarastonUVext.fits")), n.array(glob.glob( "../../LFstacks/*DEEP2*linesFittedfireflyFitsMilesUVextended.fits")) ))

stack_files.sort()

#print stack_files

for file in stack_files:
	print( file )
	if len(glob.glob(file[:-7]+"*modeled*"))>0:
		print( file, "out" )
		continue

	mm=ModelSpectraStacks(file)
	#mm.plot_fit()
	mm.fit_lines_to_lineSpectrum()
	mm.compute_derived_quantities()
	mm.save_spectrum()
	#mm.plot_spectrum()


