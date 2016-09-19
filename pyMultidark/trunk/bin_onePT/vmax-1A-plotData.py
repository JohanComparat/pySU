from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import astropy.cosmology as co
cosmo = co.Planck13

#Quantity studied
qty = "vmax"
cos = 'cen'
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits")
 )[1].data

NminCount = 1000
limits_04 = [100, 400]
limits_10 = [250, 1000]
limits_25 = [600, 1300]
limits_40 = [1200, 1600]

zmin = -0.01
zmax = 2.3

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSel = lib.nSelection(data, NminCount )
# altogether
ok = (zSel) & (mSel) & (nSel)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


# x coordinates definition
log_vmax = (n.log10(data["log_"+qty+"_min"])+n.log10(data["log_"+qty+"_max"]))/2.
vmax = 10**log_vmax
#print len(vmax), n.min(vmax), n.max(vmax)


# y coordinates
#log_VF_a = n.log10( vmax**4. * data["dNdVdlnM_"+cos][ok])
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# NOW PLOTTING ALL THE DATA
lib.plot_vmax_function_data(log_vmax[ok], log_VF[ok], data["redshift"][ok], zmin = -0.01, zmax = 2.3, cos=cos)
"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")
"""
cos = 'sat'

# y coordinates
#log_VF_a = n.log10( vmax**4. * data["dNdVdlnM_"+cos][ok])
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# NOW PLOTTING ALL THE DATA
lib.plot_vmax_function_data(log_vmax[ok], log_VF[ok], data["redshift"][ok], zmin = -0.01, zmax = 2.3, cos=cos)
"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")
"""
# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = "cen")