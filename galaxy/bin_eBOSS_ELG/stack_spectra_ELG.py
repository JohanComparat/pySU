#! /usr/bin/env python

"""
This script produces the stacks for samples defined by a list of spectra identifiers.

nohup python3 stack_spectra_ELG.py > stack_spectra_ELG.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"wwwDir/sdss/elg/stacks/KZ/v2/")

def stack_it(specList = join(spec_dir, "catalog_AGN.dat") ):
	outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
	print(specList, outfile)
	stack=sse.SpectraStackingEBOSS(specList, outfile, KZ_input=True, l_start=3.3, l_end=3.9 )
	stack.createStackMatrix()
	stack.stackSpectra()

stack_it(join(spec_dir, "catalog_lowz_liner_test.dat") )
#stack_it(join(spec_dir, "catalog_lowz_liner.dat") )
#stack_it(join(spec_dir, "catalog_lowz_AGN.dat") )
#stack_it(join(spec_dir, "catalog_lowz_comp.dat") )
#stack_it(join(spec_dir, "catalog_lowz_sf.dat") )

