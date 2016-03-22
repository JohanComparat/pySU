import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
import glob	

list_of_stacks_eb67 = glob.glob(join("/uufs/chpc.utah.edu/common/home/u0992342/eboss67/grz_stacks/","*.asc"))
list_of_stacks_eb17 = glob.glob(join("/uufs/chpc.utah.edu/common/home/u0992342/eboss17/grz_stacks/","*.asc"))

list_of_stacks = n.hstack((list_of_stacks_eb67,list_of_stacks_eb17))
list_of_stacks.sort()

for ii, el in enumerate(list_of_stacks):
		PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol = n.loadtxt(el, unpack=True)
		g_min = n.min(gmag)
		g_max = n.max(gmag)
		gr_min = n.min(grcol)
		gr_max = n.max(grcol)
		rz_min = n.min(grcol)
		rz_max = n.max(grcol)
		st=SpectraStacking(el, Nspec = 100, dLambda = 0.00005)
		outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG",el[:-4].split('/')[-1] + "_stack.fits")
		st.stackEbossPlateSpectra(PLATE.astype(int),MJD.astype(int),FIBERID.astype(int),REDSHIFT,outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)

