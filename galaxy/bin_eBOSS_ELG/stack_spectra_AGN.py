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

specList = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "xagn_0.38.asc")

print( specList )

outfile = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".stack")
stack=sse.SpectraStackingEBOSS(specList, outfile)
print(outfile)
stack.stackSpectra()

outfile = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".UVstack")
stack=sse.SpectraStackingEBOSS(specList, outfile)
print(outfile)
stack.stackSpectra_UVnormed()
	
