#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"wwwDir/sdss/elg/stacks/KZ/")

def stack_it(specList = join(spec_dir, "catalog_AGN.ascii") ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	stack=sse.SpectraStackingEBOSS(specList, outfile)
	print(outfile)
	stack.stackSpectra()

stack_it(join(spec_dir, "catalog_AGN.ascii") )
stack_it(join(spec_dir, "catalog_comp.ascii") )
stack_it(join(spec_dir, "catalog_liner.ascii") )
stack_it(join(spec_dir, "catalog_sf.ascii") )

