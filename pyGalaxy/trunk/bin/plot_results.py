#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import numpy as n
import os
from os.path import join
import astropy.io.fits as fits

G05 = n.loadtxt( join(os.environ['SPECTRASTACKS_DIR'], "biblioPoints", "gallazzi-2005.data"), unpack= True)
path_to_summary_table = join(os.environ['SPECTRASTACKS_DIR'], "results", "table_fullSpecFit_v0.VA.fits")

G05 = n.loadtxt( join(os.environ['SPECTRASTACKS_DIR'], "biblioPoints", "gallazzi-2005.data"), unpack= True)
path_to_summary_table = join(os.environ['SPECTRASTACKS_DIR'], "results", "table_fullSpecFit_v0.VA.fits")
data =fits.open(path_to_summary_table)[1].data
path_to_summary_table = join(os.environ['SPECTRASTACKS_DIR'], "results", "table_lineSpecFit_v0.VA.fits")
datL =fits.open(path_to_summary_table)[1].data


#first check that lline fits are compatible with the luminosity bin

fl = data['L_mean']/(4*n.pi*data['dL']**2.)

O2 = (data['lineWavelength']== 3728.)
O3  = (data['lineWavelength']== 5007.)
Hb = (data['lineWavelength']== 4862.) & (data['H1_4862_flux']>0)

chi2 = (data['H1_4862_flux'][Hb] - fl[Hb])/data['H1_4862_fluxErr'][Hb]

"""
ok = (data[2]<5)&(data[2]>0.)

lineWavelength,Survey,Redshift,L_MIN,L_MAX,L_MEAN,N_in_stack,R_stack,spm_light_age,spm_light_age_err_plus,spm_light_age_err_minus,spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus,spm_stellar_mass,spm_stellar_mass_err_plus,spm_stellar_mass_err_minus,spm_EBV,gp_EBV_4862_4341,gp_EBV_4862_4341_err,gp_EBV_4862_4102,gp_EBV_4862_4102_err,gp_BD_4102_4341,gp_BD_4102_4341_err,gp_SFR_O2_3728,gp_SFR_O2_3728_err,gp_SFR_H1_4862,gp_SFR_H1_4862_err,gp_12logOH_tremonti04,gp_12logOH_tremonti04_err,gp_12logOH_tremonti04_intrinsic,gp_12logOH_tremonti04_intrinsic_err = data.T[ok].T

n.savetxt(join(os.environ['SPECTRASTACKS_DIR'], "results", "table_firefly.tex"), n.transpose([Survey, lineWavelength, Redshift, N_in_stack, n.log10(n.round(L_MEAN,2)), n.round(spm_light_age,3), n.round(spm_stellar_mass,2), n.round(spm_light_metallicity,3) ]) ,delimiter = " & " ,fmt='%2.3f', newline="\\\\ ")
"""
# EBV comparison
#spm_EBV
#gp_EBV_4862_4341,gp_EBV_4862_4341_err,


p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
#ok = (datL['H1_4862_flux']>3*datL['H1_4862_fluxErr'])  & (datL['H1_4862_flux']>0) & (datL['H1_4862_fluxErr']>0)  & (datL['H1_4341_flux']>3*datL['H1_4341_fluxErr'])& (datL['H1_4341_flux']>0)&(datL['H1_4341_fluxErr']>0)& (datL['EBV_4862_4341']!=-9999.99) &  (datL['EBV_4862_4341_err']!=-9999.99)
#p.errorbar(datL['spm_EBV'][ok], datL['EBV_4862_4341'][ok], yerr=datL['EBV_4862_4341_err'][ok],fmt='o',elinewidth=2, mfc='none',label='resi')
ok = (data['H1_4862_flux']>3*data['H1_4862_fluxErr'])  & (data['H1_4862_flux']>0) & (data['H1_4862_fluxErr']>0)  & (data['H1_4341_flux']>3*data['H1_4341_fluxErr'])& (data['H1_4341_flux']>0)&(data['H1_4341_fluxErr']>0)& (data['EBV_4862_4341']!=-9999.99) &  (data['EBV_4862_4341_err']!=-9999.99)
p.errorbar(data['spm_EBV'][ok], data['EBV_4862_4341'][ok], yerr=data['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label='full')
ok = (datL['H1_4862_flux']>3*datL['H1_4862_fluxErr'])  & (datL['H1_4862_flux']>0) & (datL['H1_4862_fluxErr']>0)  & (datL['H1_4341_flux']>3*datL['H1_4341_fluxErr'])& (datL['H1_4341_flux']>0)&(datL['H1_4341_fluxErr']>0)& (datL['EBV_4862_4341']!=-9999.99) &  (datL['EBV_4862_4341_err']!=-9999.99)
p.errorbar(datL['spm_EBV'][ok], datL['EBV_4862_4341'][ok], yerr=datL['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label='GP')
p.plot([-0.1,1.5],[-0.1,1.5],'k--')
p.legend(loc=4)
p.xlabel('E(B-V) SPM')
p.ylabel(r'E(B-V) GP $H\beta -H\delta$')
p.xlim((-0.1,1.5))
p.ylim((-0.1,1.5))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-comparison-1.png"))
p.clf()



p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
ok1 = (data['H1_4862_flux']>3*data['H1_4862_fluxErr'])  & (data['H1_4862_flux']>0) & (data['H1_4862_fluxErr']>0)  & (data['H1_4341_flux']>3*data['H1_4341_fluxErr'])& (data['H1_4341_flux']>0)&(data['H1_4341_fluxErr']>0)& (data['EBV_4862_4341']!=-9999.99) &  (data['EBV_4862_4341_err']!=-9999.99)
ok2 = (data['H1_4862_flux']>3*data['H1_4862_fluxErr'])  & (data['H1_4862_flux']>0) & (data['H1_4862_fluxErr']>0)  & (data['H1_4102_flux']>3*data['H1_4102_fluxErr'])& (data['H1_4102_flux']>0)&(data['H1_4102_fluxErr']>0)& (data['EBV_4862_4102']!=-9999.99) &  (data['EBV_4862_4102_err']!=-9999.99)
ok=(ok1)&(ok2)
p.errorbar(data['EBV_4862_4102'][ok], data['EBV_4862_4341'][ok], xerr=data['EBV_4862_4102_err'][ok], yerr=data['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label='full')
ok1 = (datL['H1_4862_flux']>3*datL['H1_4862_fluxErr'])  & (datL['H1_4862_flux']>0) & (datL['H1_4862_fluxErr']>0)  & (datL['H1_4341_flux']>3*datL['H1_4341_fluxErr'])& (datL['H1_4341_flux']>0)&(datL['H1_4341_fluxErr']>0)& (datL['EBV_4862_4341']!=-9999.99) &  (datL['EBV_4862_4341_err']!=-9999.99)
ok2 = (datL['H1_4862_flux']>3*datL['H1_4862_fluxErr'])  & (datL['H1_4862_flux']>0) & (datL['H1_4862_fluxErr']>0)  & (datL['H1_4102_flux']>3*datL['H1_4102_fluxErr'])& (datL['H1_4102_flux']>0)&(datL['H1_4102_fluxErr']>0)& (datL['EBV_4862_4102']!=-9999.99) &  (datL['EBV_4862_4102_err']!=-9999.99)
ok=(ok1)&(ok2)
p.errorbar(datL['EBV_4862_4102'][ok], datL['EBV_4862_4341'][ok], xerr=datL['EBV_4862_4102_err'][ok], yerr=datL['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label='GP')
p.plot([-0.1,1.5],[-0.1,1.5],'k--')
p.legend(loc=4)
p.xlabel(r'E(B-V) GP $H\beta -H\gamma$')
p.ylabel(r'E(B-V) GP $H\beta -H\delta$')
p.xlim((-0.1,1.5))
p.ylim((-0.1,1.5))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-comparison-2.png"))
p.clf()


Nst = (data['N_in_stack']==200)&(data['Survey']==2)
EBvok = (Nst)&(data['H1_4862_flux']>data['H1_4862_fluxErr'])  & (data['H1_4862_flux']>0) & (data['H1_4862_fluxErr']>0)  & (data['H1_4341_flux']>data['H1_4341_fluxErr'])& (data['H1_4341_flux']>0)&(data['H1_4341_fluxErr']>0)& (data['EBV_4862_4341']!=-9999.99) &  (data['EBV_4862_4341_err']!=-9999.99)
O2 = (data['lineWavelength']== 3728.) & (EBvok)
O3  = (data['lineWavelength']== 5007.)& (EBvok)
Hb = (data['lineWavelength']== 4862.) & (EBvok)
data['Redshift'][O2]


p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
for zz in set(data['Redshift'][O2]):
	ok = (O2)&(data['Redshift']==zz)
	p.errorbar(data['L_mean'][ok], data['EBV_4862_4341'][ok], xerr=[data['L_mean'][ok]-data['L_min'][ok],-data['L_mean'][ok]+data['L_max'][ok]],yerr=data['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label=str(n.round(zz,3)))

p.legend(loc=0)
p.xlabel('[OII] line Luminosity')
p.ylabel(r'E(B-V) GP $H\beta -H\delta$')
p.xscale('log')
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-line-o2.png"))
p.clf()


Nst = (datL['N_in_stack']==200)&(datL['Survey']==2)
EBvok = (Nst)&(datL['H1_4862_flux']>datL['H1_4862_fluxErr'])  & (datL['H1_4862_flux']>0) & (datL['H1_4862_fluxErr']>0)  & (datL['H1_4341_flux']>datL['H1_4341_fluxErr'])& (datL['H1_4341_flux']>0)&(datL['H1_4341_fluxErr']>0)& (datL['EBV_4862_4341']!=-9999.99) &  (datL['EBV_4862_4341_err']!=-9999.99)
O2 = (datL['lineWavelength']== 3728.) & (EBvok)
O3  = (datL['lineWavelength']== 5007.)& (EBvok)
Hb = (datL['lineWavelength']== 4862.) & (EBvok)
datL['Redshift'][O2]


p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
for zz in set(datL['Redshift'][O2]):
	ok = (O2)&(datL['Redshift']==zz)
	p.errorbar(datL['L_mean'][ok], datL['EBV_4862_4341'][ok], xerr=[datL['L_mean'][ok]-datL['L_min'][ok],-datL['L_mean'][ok]+datL['L_max'][ok]],yerr=datL['EBV_4862_4341_err'][ok],fmt='o',elinewidth=1, mfc='none',label=str(n.round(zz,3)))

p.legend(loc=0)
p.xlabel('[OII] line Luminosity')
p.ylabel(r'E(B-V) GP $H\beta -H\delta$')
p.xscale('log')
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-line-o2-gp.png"))
p.clf()


















p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
ok = (data['logOH_tremonti04']!=-9999.99) &  (data['logOH_tremonti04_err']!=-9999.99)
p.errorbar(x=data['spm_light_metallicity'][ok], xerr=[data['spm_light_metallicity_err_plus'][ok], data['spm_light_metallicity_err_minus'][ok]], y=data['logOH_tremonti04'][ok], yerr=data['logOH_tremonti04_err'][ok],fmt='o',elinewidth=2, mfc='none',label='full')
ok = (datL['logOH_tremonti04']!=-9999.99) &  (datL['logOH_tremonti04_err']!=-9999.99)
p.errorbar(x=datL['spm_light_metallicity'][ok], xerr=[datL['spm_light_metallicity_err_plus'][ok], datL['spm_light_metallicity_err_minus'][ok]], y=datL['logOH_tremonti04'][ok], yerr=datL['logOH_tremonti04_err'][ok],fmt='o',elinewidth=2, mfc='none',label='GP')
#p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('log(Z/Zsun) SPM')
p.ylabel(r'12+log(O/H) (line ratios Tremonti 04 estimator)')
#p.xlim((-1,2))
#p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "metal-comparison-1.png"))
p.clf()

# SFR comparison
o2_sfr_conv =  10**(0.27) * 10**(-41) 
hb_sfr_conv = 10**(0.58) * 10**(-41)

sfrO2 = data['L_mean']*o2_sfr_conv
sfrO2_up = data['L_max']*o2_sfr_conv
sfrO2_low = data['L_min']*o2_sfr_conv

sfrHb = data['L_mean']*hb_sfr_conv
sfrHb_up = data['L_max']*hb_sfr_conv
sfrHb_low = data['L_min']*hb_sfr_conv


Nst = (data['N_in_stack']==400)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=data['Redshift'][hb],label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=data['Redshift'][o2],label='O2')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o3], data['spm_light_metallicity'][o3], marker ='^',  s=30, c=data['Redshift'][o3],label='O3')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('redshift')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-redshift-all-hb-o2-o3-400.png"))
p.clf()


Nst = (data['N_in_stack']==100)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=data['Redshift'][hb],label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=data['Redshift'][o2],label='O2')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o3], data['spm_light_metallicity'][o3], marker ='^',  s=30, c=data['Redshift'][o3],label='O3')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('redshift')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-redshift-all-hb-o2-o3-100.png"))
p.clf()


Nst = (data['N_in_stack']==200)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=data['Redshift'][hb],label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=data['Redshift'][o2],label='O2')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o3], data['spm_light_metallicity'][o3], marker ='^',  s=30, c=data['Redshift'][o3],label='O3')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('redshift')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-redshift-all-hb-o2-o3-200.png"))
p.clf()



Nst = (data['N_in_stack']==200)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-all-hb-o2-200.png"))
p.clf()


Nst = (data['N_in_stack']==400)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-all-hb-o2-400.png"))
p.clf()


Nst = (data['N_in_stack']==100)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-all-hb-o2-100.png"))
p.clf()






Nst = (data['N_in_stack']==100)&(data['Redshift']>=0.7)&(data['Redshift']<=0.85)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-07z085-hb-o2-100.png"))
p.clf()


Nst = (data['N_in_stack']==200)&(data['Redshift']>=0.7)&(data['Redshift']<=0.85)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-07z085-hb-o2-200.png"))
p.clf()



Nst = (data['N_in_stack']==400)&(data['Redshift']>=0.7)&(data['Redshift']<=0.85)
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ttt = p.scatter(data['spm_stellar_mass'][hb], data['spm_light_metallicity'][hb], marker='s', s=30, c=n.log10(sfrHb[hb]),label='Hb')
ttt.set_edgecolor('face')
ttt = p.scatter(data['spm_stellar_mass'][o2], data['spm_light_metallicity'][o2], marker ='o',  s=30, c=n.log10(sfrO2[o2]),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "mass-metallicity-sfr-07z085-hb-o2-400.png"))
p.clf()




# mass metallicity relation
# Mannucci et al. 2010 :
y_rel = lambda x_rel : 8.90 + 0.39*x_rel - 0.20*x_rel*x_rel - 0.077*x_rel*x_rel*x_rel + 0.064*x_rel*x_rel*x_rel*x_rel
x_rel = n.arange(8.5,11.5,0.05) - 10 
y_pr = y_rel(x_rel)

Nst =  (data['logOH_tremonti04']!=-9999.99) # (data['Redshift']>=0.7)&(data['Redshift']<=0.85) &(data['N_in_stack']==400)&
hb = (data['lineWavelength']==4862.0) & (Nst)
o2 = (data['lineWavelength']==3728.0) & (Nst)
o3 = (data['lineWavelength']==5007.0) & (Nst)

x_obs_1 = data['spm_stellar_mass'] - 0.32 * n.log10(sfrHb) - 10
x_obs_2 = data['spm_stellar_mass'] - 0.32 * n.log10(sfrO2) - 10
y_obs_1 = 9.185-0.313* data['logOH_tremonti04'] - 0.264 * data['logOH_tremonti04']**2 - 0.321 * data['logOH_tremonti04']**3

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
p.plot(x_obs_1[hb], data['logOH_tremonti04'][hb], 'bo',label='Hb')
p.plot(x_obs_2[o2], data['logOH_tremonti04'][o2], 'r*',label='O2')
p.plot(x_rel, y_pr, label='Manucci 2010')
p.ylabel(r'$12+log(OH)$')
p.xlabel('$\log(M_*)-0.32 log(SFR) - 10$')
p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=2)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "manucci-2010-relation-all-hb-o2.png"))
p.clf()


Nst =  (datL['logOH_tremonti04']!=-9999.99) # (datL['Redshift']>=0.7)&(datL['Redshift']<=0.85) &(datL['N_in_stack']==400)&
hb = (datL['lineWavelength']==4862.0) & (Nst)
o2 = (datL['lineWavelength']==3728.0) & (Nst)
o3 = (datL['lineWavelength']==5007.0) & (Nst)

x_obs_1 = datL['spm_stellar_mass'] - 0.32 * n.log10(sfrHb) - 10
x_obs_2 = datL['spm_stellar_mass'] - 0.32 * n.log10(sfrO2) - 10
y_obs_1 = 9.185-0.313* datL['logOH_tremonti04'] - 0.264 * datL['logOH_tremonti04']**2 - 0.321 * datL['logOH_tremonti04']**3

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
p.plot(x_obs_1[hb], datL['logOH_tremonti04'][hb], 'bo',label='Hb')
p.plot(x_obs_2[o2], datL['logOH_tremonti04'][o2], 'r*',label='O2')
p.plot(x_rel, y_pr, label='Manucci 2010')
p.ylabel(r'$12+log(OH)$')
p.xlabel('$\log(M_*)-0.32 log(SFR) - 10$')
p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=2)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "manucci-2010-GP-relation-all-hb-o2.png"))
p.clf()




















p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
p.errorbar(x=gp_SFR_O2_3728[ok], xerr = gp_SFR_O2_3728_err[ok], y=gp_SFR_H1_4862[ok], yerr=gp_SFR_H1_4862_err[ok],fmt='o',elinewidth=2, mfc='none')
p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('SFR GP [OII]')
p.ylabel(r'SFR GP H$\beta$')
#p.xlim((-1,2))
#p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "sfr-comparison-2.png"))
p.clf()

# mass metallicity relation
# Mannucci et al. 2010 :
y_rel = lambda x_rel : 8.90 + 0.39*x_rel - 0.20*x_rel*x_rel - 0.077*x_rel*x_rel*x_rel + 0.064*x_rel*x_rel*x_rel*x_rel

x_obs_1 = spm_stellar_mass - 0.32 * n.log10(gp_SFR_O2_3728) - 10
x_obs_2 = spm_stellar_mass - 0.32 * n.log10(gp_SFR_H1_4862) - 10
y_obs_1 = 9.185-0.313*spm_light_metallicity - 0.264 *spm_light_metallicity**2 - 0.321 *spm_light_metallicity**3

x_rel = n.arange(8.5,11.5,0.05) - 10 
y_pr = y_rel(x_rel)

#base = (gp_SFR_O2_3728>0)& (gp_SFR_O2_3728_err>0)&(spm_stellar_mass >0)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
p.plot(G05[0],G05[1], 'k',label ="G05")
p.plot(G05[0],G05[2], 'k--')
p.plot(G05[0],G05[3], 'k--')

#p.plot(x_obs_1, y_obs_1, )
ok = (lineWavelength==4862.)
ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker='s', s=30, c=n.log10(L_MEAN[ok]*10**(0.58) * 10**(-41) ),label='Hb')
ttt.set_edgecolor('face')
#ok = (base) &(lineWavelength==5007.)
#ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker ='o',  s=30, c=n.log10(gp_SFR_H1_4862[ok]),label='O3')
#ttt.set_edgecolor('face')
ok = (lineWavelength==3728.)
ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker = '*', s=30, c=n.log10(L_MEAN[ok]*10**(0.27) * 10**(-41) ),label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('log(SFR)')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "sfr-mass-z-allredshift.png"))
p.clf()


p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.75])
#p.plot(x_obs_1, y_obs_1, )
ok = (lineWavelength==4862.)
ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker='s', s=30, c=Redshift[ok],label='Hb')
ttt.set_edgecolor('face')
ok = (lineWavelength==5007.)
ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker ='o',  s=30, c=Redshift[ok],label='O3')
ttt.set_edgecolor('face')
ok = (lineWavelength==3728.)
ttt = p.scatter(spm_stellar_mass[ok], spm_light_metallicity[ok], marker = '*', s=30, c=Redshift[ok],label='O2')
ttt.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label('redshift')
p.ylabel(r'$log([Z/H])$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "sfr-mass-redshift-allredshift.png"))
p.clf()