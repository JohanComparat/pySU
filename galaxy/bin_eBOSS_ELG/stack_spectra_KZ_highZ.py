#! /usr/bin/env python

"""
This script produces the stacks for samples defined by a list of spectra identifiers.

nohup python2 stack_spectra_KZ_highZ.py > stack_spectra_KZ_highZ.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"wwwDir/sdss/elg/stacks/KZ/high_z/")

def stack_it(specList = join(spec_dir, "catalog_AGN.ascii") ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	print(specList, outfile)
	stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
	stack.createStackMatrix()
	stack.stackSpectra()

stack_it(join(spec_dir, 'catalog_AGN.ascii'   ) )
stack_it(join(spec_dir, 'catalog_liner.ascii' ) )
stack_it(join(spec_dir, 'catalog_comp.ascii'  ) )
stack_it(join(spec_dir, 'catalog_sf.ascii'    ) )





