import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
import glob	

dir_to_stack_list = ""
list_of_stacks = glob.glob(join(dir_to_stack_list,"stackList*.fits"))

list_of_stacks = n.array(["/uufs/chpc.utah.edu/common/home/u0936736/stack_SDSSMAIN/specIDS_z_wfc_wcomp_SDSS_test.fits"])

el = list_of_stacks[0]
st=SpectraStacking(el, Nspec = 998, dLambda = 0.00005)
outPutFileName = el[:-5] + "_stack.fits"
st.stackSdssMainSpectra(outPutFileName)

"""
for ii, el in enumerate(list_of_stacks):
		st=SpectraStacking(el, Nspec = 100, dLambda = 0.00005)
		outPutFileName = el[:-5] + "_stack.fits"
		st.stackSdssMainSpectra(self,outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)

"""
		
