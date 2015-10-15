from InterpretStacksAndLFs import *

lf_files=n.hstack(( n.array(glob.glob( "/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/??_????/*-z?.???.fits")), n.array(glob.glob( "/home/comparat/database/VVDS/products/emissionLineLuminosityFunctions/??_????/*-z?.???.fits")) ))

lf_files.sort()

for file in lf_files :
	print file
	mm = IntrinsicStacksAndLFs( file )
	mm.apply_correction()
