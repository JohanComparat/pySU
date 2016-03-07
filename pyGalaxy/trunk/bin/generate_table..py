
from ModelSpectraStacks import *
import glob
from os.path import join
import os

stack_files=n.array(glob.glob( join( os.environ['SPECTRASTACKS_DIR'], "data", "??_????", "*DEEP2*.fits" ) ) )
stack_files.sort()

#print stack_files

file = stack_files[0]
for file in stack_files:
	mm=ModelSpectraStacks(file)
	print file.split('/')[-1], n.round(mm.wl.min()), " & " , n.round(mm.wl.max())