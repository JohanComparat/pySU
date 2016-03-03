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

path_to_summary_table = join(os.environ['SPECTRASTACKS_DIR'], "results", "table_v0.data")

data = n.loadtxt(path_to_summary_table)	

ok = (data[2]<5)&(data[2]>0.)

lineWavelength,Survey,Redshift,L_MIN,L_MAX,L_MEAN,N_in_stack,R_stack,spm_light_age,spm_light_age_err_plus,spm_light_age_err_minus,spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus,spm_stellar_mass,spm_stellar_mass_err_plus,spm_stellar_mass_err_minus,spm_EBV,gp_EBV_4862_4341,gp_EBV_4862_4341_err,gp_EBV_4862_4102,gp_EBV_4862_4102_err,gp_BD_4102_4341,gp_BD_4102_4341_err,gp_SFR_O2_3728,gp_SFR_O2_3728_err,gp_SFR_H1_4862,gp_SFR_H1_4862_err,gp_12logOH_tremonti04,gp_12logOH_tremonti04_err,gp_12logOH_tremonti04_intrinsic,gp_12logOH_tremonti04_intrinsic_err = data.T[ok].T

n.savetxt(join(os.environ['SPECTRASTACKS_DIR'], "results", "table_firefly.tex"), n.transpose([Survey, lineWavelength, Redshift, N_in_stack, n.log10(n.round(L_MEAN,2)), n.round(spm_light_age,3), n.round(spm_stellar_mass,2), n.round(spm_light_metallicity,3) ]) ,delimiter = " & " ,fmt='%2.3f', newline="\\\\ ")

# EBV comparison
#spm_EBV
#gp_EBV_4862_4341,gp_EBV_4862_4341_err,

ok = (gp_EBV_4862_4341!=-9999.99) &  (gp_EBV_4862_4341_err!=-9999.99)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
p.errorbar(spm_EBV[ok], gp_EBV_4862_4341[ok], yerr=gp_EBV_4862_4341_err[ok],fmt='o',elinewidth=2, mfc='none')
p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('E(B-V) SPM')
p.ylabel(r'E(B-V) GP $H\beta -H\delta$')
p.xlim((-1,2))
p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-comparison-1.png"))
p.clf()

#spm_EBV
#gp_EBV_4862_4102,gp_EBV_4862_4102_err

ok = (gp_EBV_4862_4102!=-9999.99) &  (gp_EBV_4862_4102_err!=-9999.99)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
p.errorbar(spm_EBV[ok], gp_EBV_4862_4102[ok], yerr=gp_EBV_4862_4102_err[ok],fmt='o',elinewidth=2, mfc='none')
p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('E(B-V) SPM')
p.ylabel(r'E(B-V) GP $H\beta -H\gamma$')
p.xlim((-1,2))
p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "ebv-comparison-2.png"))
p.clf()

# metallicity comparison
#spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus
#gp_12logOH_tremonti04,gp_12logOH_tremonti04_err,

ok = (gp_12logOH_tremonti04!=-9999.99) &  (gp_12logOH_tremonti04_err!=-9999.99)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
p.errorbar(x=spm_light_metallicity[ok], xerr=[spm_light_metallicity_err_plus[ok], spm_light_metallicity_err_minus[ok]], y=gp_12logOH_tremonti04[ok], yerr=gp_12logOH_tremonti04_err[ok],fmt='o',elinewidth=2, mfc='none')
p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('log(Z/Zsun) SPM')
p.ylabel(r'12+log(O/H) GP')
#p.xlim((-1,2))
#p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "metal-comparison-1.png"))
p.clf()

#spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus
#gp_12logOH_tremonti04_intrinsic,gp_12logOH_tremonti04_intrinsic_err

ok = (gp_12logOH_tremonti04_intrinsic!=-9999.99) &  (gp_12logOH_tremonti04_intrinsic_err!=-9999.99)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.8,0.8])
p.errorbar(x=spm_light_metallicity[ok], xerr=[spm_light_metallicity_err_plus[ok], spm_light_metallicity_err_minus[ok]], y=gp_12logOH_tremonti04_intrinsic[ok], yerr=gp_12logOH_tremonti04_intrinsic_err[ok],fmt='o',elinewidth=2, mfc='none')
p.plot([-1,2],[-1,2],'k--')
p.legend(loc=2)
p.xlabel('log(Z/Zsun) SPM')
p.ylabel(r'12+log(O/H) GP intrinsic')
#p.xlim((-1,2))
#p.ylim((-1,2))
p.grid()
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "metal-comparison-2.png"))
p.clf()

# SFR comparison

ok = (gp_SFR_O2_3728!=-9999.99) &  (gp_SFR_O2_3728_err!=-9999.99)&(gp_SFR_H1_4862!=-9999.99) &  (gp_SFR_H1_4862_err!=-9999.99)

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
p.ylabel(r'$[Z/H]$')
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
p.ylabel(r'$[Z/H]$')
p.xlabel('$\log(M_*)$')
#p.title(r'$0.7<z<0.85$ ELG')
p.legend(loc=3)
p.savefig( join(os.environ['SPECTRASTACKS_DIR'], "plots", "sfr-mass-redshift-allredshift.png"))
p.clf()