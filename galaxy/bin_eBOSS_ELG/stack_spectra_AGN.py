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


spec_dir = join(os.environ['HOME'],"wwwDir/sdss/agn/stacks")

def stack_it(specList = join(spec_dir, "catalog.ascii") ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	stack=sse.SpectraStackingEBOSS(specList, outfile)
	stack.createStackMatrix()
	print(outfile)
	stack.stackSpectra()

stack_it(join(spec_dir, "s82xagn-gal_ell.asc") )
stack_it(join(spec_dir, "s82xagn-gal_SB.asc") )
stack_it(join(spec_dir, "s82xagn-gal_spi.asc") )
stack_it(join(spec_dir, "s82xagn-stars.asc") )
stack_it(join(spec_dir, "s82xagn-t1.asc") )
stack_it(join(spec_dir, "s82xagn-t2.asc") )

