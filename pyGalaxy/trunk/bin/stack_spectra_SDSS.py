import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
import glob	

dir_to_stack_list = ""
list_of_stacks = glob.glob(join(dir_to_stack_list,"stackList*.fits"))

for ii, el in enumerate(list_of_stacks):
		g_min = 
		g_max =
		gr_min = 
		gr_max = 
		rz_min = 
		rz_max = 
		st=SpectraStacking(el, Nspec = 100, dLambda = 0.00005)
		outPutFileName = el[:-5] + "_stack.fits"
		st.stackSdssSpectra(outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)

