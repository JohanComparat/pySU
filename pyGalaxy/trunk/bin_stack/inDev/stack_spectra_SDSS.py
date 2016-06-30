import sys
import os 
from os.path import join
from SpectraStackingSDSSOnly import *
import glob	

hdus_eb67 = fits.open("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG/elg270_eboss67summaryTable_stack_comparison.fits")
hdus_eb17 = fits.open("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG/elg270_eboss17summaryTable_stack_comparison.fits")
hdus_ELG = fits.open("/uufs/chpc.utah.edu/common/home/u0936736/ELG_catalogs/spAll-ELG-v6.5.4.fits")

#'PLATE', 'MJD', 'FIBER', 'gmag', 'rzcol', 'grcol', 'Z_1', 'Z_2', 'Z_3', 'Z_ERR_1', 'Z_ERR_2', 'Z_ERR_3', 'RCHI2_1', 'RCHI2_2', 'RCHI2_3', 'CLASS_1', 'CLASS_2', 'CLASS_3'

ZQgrid  = [-2.5, -1.5, -0.5, 0.5, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75]
ZCgrid = [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
Zgrid = n.hstack(( n.arange(0,0.6,0.1), n.arange(0.6, 1.1, 0.05), n.arange(1.1, 1.7, 0.1) ))

def produce_stacks_zQ_zCont_Z(table, ZQgrid, ZCgrid, Zgrid, nameRoot="elg270_eboss67_"):
	print table.dtype
	index_g = n.ones_like(table['gmag'])*-1
	index_rz = n.ones_like(table['gmag'])*-1
	index_gr = n.ones_like(table['gmag'])*-1
	for i in range(len(ZQgrid)-1):
		for j in range(len(ZCgrid)-1):
			for k in range(len(Zgrid)-1):
				sel = (table['zQ']>ZQgrid[i])&(table['zQ']<ZQgrid[i+1]) & (table['zCont']>ZCgrid[j])&(table['zCont']<ZCgrid[j+1]) & (table['Z']>Zgrid[k])&(table['Z']<Zgrid[k+1])# &(table['Z_1']>0)&(table['Z_1']>table['Z_ERR_1'])&(table['Z_ERR_1']>0)
				index_g[sel] = i*n.ones_like(index_g[sel])
				index_rz[sel] = j*n.ones_like(index_g[sel])
				index_gr[sel] = k*n.ones_like(index_g[sel])
				
				PLATE ,   MJD  ,  FIBERID ,   REDSHIFT , z, zq, zc = table['PLATE'][sel], table['MJD'][sel], table['FIBER'][sel], table['Z'][sel], table['zQ'][sel], table['zCont'][sel]
				zq_min = n.min(zq)
				zq_max = n.max(zq)
				zc_min = n.min(zc)
				zc_max = n.max(zc)
				z_min = n.min(z)
				z_max = n.max(z)
				st=SpectraStacking("-", Nspec = 100, dLambda = 0.00005)
				suffix = "_zQ_"+str(n.round(ZQgrid[i],1))+"_zC_"+str(n.round(ZCgrid[j],1))+"_Z_"+str(n.round(Zgrid[k],1))
				outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + suffix + "_stack.fits")
				st.stackEbossPlateSpectra(PLATE.astype(int),MJD.astype(int),FIBERID.astype(int),REDSHIFT,outPutFileName, g_min = zq_min,g_max=zq_max, gr_min=zc_min, gr_max=zc_max, rz_min= z_min, rz_max = z_max)

	summaryTableName =join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + "summaryTable_stack.fits")
	col_index_g = fits.Column(name="index_g",format="I", array= index_g.astype(int))
	col_index_rz = fits.Column(name="index_rz",format="I", array= index_rz.astype(int))
	col_index_gr = fits.Column(name="index_gr",format="I", array= index_gr.astype(int))
	cols = table.columns + col_index_g + col_index_rz + col_index_gr
	tbhdu = fits.BinTableHDU.from_columns(cols)
	prihdr = fits.Header()
	prihdr['chunk'] = nameRoot
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	os.system('rm '+summaryTableName)
	thdulist.writeto(summaryTableName)


produce_stacks_zQ_zCont_Z(hdus_ELG[1].data, nameRoot="elg_zQzCz")


sys.exit()

def produce_stacks_z(table, nameRoot="elg270_eboss67"):
	print table.dtype
	zarr = table['Z_1'][((table['Z_1']>0)&(table['Z_1']>table['Z_ERR_1'])&(table['Z_ERR_1']>0))&((table['CLASS_1']=="QSO")|(table['CLASS_1']=="GALAXY"))]
	zarr.sort()
	grid  = zarr[::50]
	index_Z1 = n.ones_like(table['gmag'])*-1
	for i in range(len(grid)-1):
		sel = ((table['Z_1']>=grid[i])&(table['Z_1']<grid[i+1])&(table['Z_1']>0)&(table['Z_1']>table['Z_ERR_1'])&(table['Z_ERR_1']>0))&((table['CLASS_1']=="QSO")|(table['CLASS_1']=="GALAXY"))
		index_Z1[sel] = i*n.ones_like(index_Z1[sel])
		PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol = table['PLATE'][sel], table['MJD'][sel], table['FIBER'][sel], table['Z_1'][sel], table['gmag'][sel], table['rzcol'][sel], table['grcol'][sel]
		g_min = n.min(gmag)
		g_max = n.max(gmag)
		gr_min = n.min(grcol)
		gr_max = n.max(grcol)
		rz_min = n.min(grcol)
		rz_max = n.max(grcol)
		st=SpectraStacking("-", Nspec = 100, dLambda = 0.00005)
		suffix = "_Z1_"+str(n.round(grid[i],3))+"_"+str(n.round(grid[i+1],3))
		outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + suffix + "_stackN50.fits")
		st.stackEbossPlateSpectra(PLATE.astype(int),MJD.astype(int),FIBERID.astype(int),REDSHIFT,outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)

	summaryTableName =join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + "_summaryTable_Zstack.fits")
	col_index_Z1 = fits.Column(name="index_Z1",format="I", array= index_Z1.astype(int))
	cols = table.columns + col_index_Z1
	tbhdu = fits.BinTableHDU.from_columns(cols)
	prihdr = fits.Header()
	prihdr['chunk'] = nameRoot
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	os.system('rm '+summaryTableName)
	thdulist.writeto(summaryTableName)


produce_stacks_z(hdus_eb17[1].data, nameRoot="elg270_eboss17")
produce_stacks_z(hdus_eb67[1].data, nameRoot="elg270_eboss67")

sys.exit()

ggrid  = [21.8,22.5,22.8]
rzgrid = [0.0,0.8,1.0,2.0]
grgrid = [0.0,0.4,0.6,1.0]

def produce_stacks(table, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss67_"):
	print table.dtype
	index_g = n.ones_like(table['gmag'])*-1
	index_rz = n.ones_like(table['gmag'])*-1
	index_gr = n.ones_like(table['gmag'])*-1
	for i in range(len(ggrid)-1):
		for j in range(len(rzgrid)-1):
			for k in range(len(grgrid)-1):
				sel = (table['gmag']>ggrid[i])&(table['gmag']<ggrid[i+1]) & (table['rzcol']>rzgrid[j])&(table['rzcol']<rzgrid[j+1]) & (table['grcol']>grgrid[k])&(table['grcol']<grgrid[k+1])&(table['Z_1']>0)&(table['Z_1']>table['Z_ERR_1'])&(table['Z_ERR_1']>0)
				index_g[sel] = i*n.ones_like(index_g[sel])
				index_rz[sel] = j*n.ones_like(index_g[sel])
				index_gr[sel] = k*n.ones_like(index_g[sel])
				
				PLATE ,   MJD  ,  FIBERID ,   REDSHIFT   , gmag ,   rzcol  ,  grcol = table['PLATE'][sel], table['MJD'][sel], table['FIBER'][sel], table['Z_1'][sel], table['gmag'][sel], table['rzcol'][sel], table['grcol'][sel]
				g_min = n.min(gmag)
				g_max = n.max(gmag)
				gr_min = n.min(grcol)
				gr_max = n.max(grcol)
				rz_min = n.min(grcol)
				rz_max = n.max(grcol)
				st=SpectraStacking("-", Nspec = 100, dLambda = 0.00005)
				suffix = "_g_"+str(n.round(ggrid[i],1))+"_rz_"+str(n.round(rzgrid[j],1))+"_gr_"+str(n.round(grgrid[k],1))
				outPutFileName = join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + suffix + "_stack.fits")
				st.stackEbossPlateSpectra(PLATE.astype(int),MJD.astype(int),FIBERID.astype(int),REDSHIFT,outPutFileName, g_min = g_min,g_max=g_max, gr_min=gr_min, gr_max=gr_max, rz_min= rz_min, rz_max = rz_max)

	summaryTableName =join("/uufs/chpc.utah.edu/common/home/u0936736/stack_eBOSSELG", nameRoot + "summaryTable_stack.fits")
	col_index_g = fits.Column(name="index_g",format="I", array= index_g.astype(int))
	col_index_rz = fits.Column(name="index_rz",format="I", array= index_rz.astype(int))
	col_index_gr = fits.Column(name="index_gr",format="I", array= index_gr.astype(int))
	cols = table.columns + col_index_g + col_index_rz + col_index_gr
	tbhdu = fits.BinTableHDU.from_columns(cols)
	prihdr = fits.Header()
	prihdr['chunk'] = nameRoot
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	os.system('rm '+summaryTableName)
	thdulist.writeto(summaryTableName)


produce_stacks(hdus_eb17[1].data, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss17")
produce_stacks(hdus_eb67[1].data, ggrid, rzgrid, grgrid, nameRoot="elg270_eboss67")
