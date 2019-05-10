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

#spec_dir = join(os.environ['HOME'],"wwwDir/sdss/elg/stacks/KZ/high_z/")
spec_dir = join(os.environ['HOME'],"wwwDir/sdss/elg/stacks/KZ/low_z/")

def stack_it(specList = join(spec_dir, "catalog_AGN.ascii") ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	print(specList, outfile)
	stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
	stack.createStackMatrix()
	stack.stackSpectra()


stack_it( join(spec_dir, 'catalog_lowz_AGN.dat')   )
stack_it( join(spec_dir, 'catalog_lowz_comp.dat')  )
stack_it( join(spec_dir, 'catalog_lowz_liner.dat') )
stack_it( join(spec_dir, 'catalog_lowz_sf.dat')    )

sys.exit()

stack_it(join(spec_dir, 'catalog_AGN.ascii'   ) )
stack_it(join(spec_dir, 'catalog_liner.ascii' ) )
stack_it(join(spec_dir, 'catalog_comp.ascii'  ) )
stack_it(join(spec_dir, 'catalog_sf.ascii'    ) )

spec_dir = join(os.environ['HOME'],"afs_comparat/sdss/elg/stacks/KZ/high_z/")

specList = join(spec_dir, 'catalog_AGN.ascii'   )
outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
print(specList, outfile)
stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack.stackSpectra()

specList = join(spec_dir, 'catalog_liner.ascii'   )
outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
print(specList, outfile)
stack2=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack2.stackSpectra()

specList = join(spec_dir, 'catalog_comp.ascii'   )
outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
print(specList, outfile)
stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack.stackSpectra()

specList = join(spec_dir, 'catalog_sf.ascii'   )
outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
print(specList, outfile)
stack2=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack2.stackSpectra()

spec_dir = join(os.environ['HOME'],"afs_comparat/sdss/elg/stacks/KZ/low_z/")

specList = join(spec_dir, 'catalog_lowz_AGN.dat')                         
outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
print(specList, outfile)
stack=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack.stackSpectra()

specList = join(spec_dir, 'catalog_lowz_comp.dat')                         
outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
print(specList, outfile)
stack2=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack2.stackSpectra()

specList = join(spec_dir, 'catalog_lowz_liner.dat')                         
outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
print(specList, outfile)
stack3=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack3.stackSpectra()

specList = join(spec_dir, 'catalog_lowz_sf.dat')
outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
print(specList, outfile)
stack4=sse.SpectraStackingEBOSS(specList, outfile, l_start=3.3, l_end=3.9 )
stack4.stackSpectra()
