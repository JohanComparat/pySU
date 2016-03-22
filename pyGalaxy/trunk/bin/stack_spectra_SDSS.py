import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
import glob	

hdus_eb67 = fits.open("/uufs/chpc.utah.edu/common/home/u0992342/eboss67/elg270_eboss67_3zbest.fits")
hdus_eb17 = fits.open("/uufs/chpc.utah.edu/common/home/u0992342/eboss17/elg270_eboss17_3zbest.fits")

#'PLATE', 'MJD', 'FIBER', 'gmag', 'rzcol', 'grcol', 'Z_1', 'Z_2', 'Z_3', 'Z_ERR_1', 'Z_ERR_2', 'Z_ERR_3', 'RCHI2_1', 'RCHI2_2', 'RCHI2_3', 'CLASS_1', 'CLASS_2', 'CLASS_3'

ggrid  = [21.8,22.5,22.8]
rzgrid = [0.0,0.8,1.0,2.0]
grgrid = [0.0,0.4,0.6,1.0]

def produce_stacks(table, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss67_"):
	print table.dtype
	for i in range(len(ggrid)-1):
		for j in range(len(rzgrid)-1):
			for k in range(len(grgrid)-1):
				sel = (table['gmag']>ggrid[i])&(table['gmag']<ggrid[i+1]) & (table['rzcol']>rzgrid[j])&(table['rzcol']<rzgrid[j+1]) & (table['grcol']>grgrid[k])&(table['grcol']<grgrid[k+1])&(table['Z_1']>0)&(table['Z_1']>table['Z_ERR_1'])&(table['Z_ERR_1']>0)
				PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol = table['PLATE'][sel], table['MJD'][sel], table['FIBER'][sel], table['Z_1'][sel], table['gmag'][sel], table['rzcol'][sel], table['grcol'][sel]
				g_min = n.min(gmag)
				g_max = n.max(gmag)
				gr_min = n.min(grcol)
				gr_max = n.max(grcol)
				rz_min = n.min(grcol)
				rz_max = n.max(grcol)
				st=SpectraStacking("-", Nspec = 100, dLambda = 0.00005)
				suffix = "_g_"+str(n.round(ggrid[i],1))+"_rz_"+str(n.round(rzgrid[i],1))+"_gr_"+str(n.round(grgrid[i],1))
				outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + suffix + "_stack.fits")
				st.stackEbossPlateSpectra(PLATE.astype(int),MJD.astype(int),FIBERID.astype(int),REDSHIFT,outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)


produce_stacks(hdus_eb17[1].data, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss17_")
produce_stacks(hdus_eb67[1].data, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss67_")
