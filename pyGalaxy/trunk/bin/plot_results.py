#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

import os
from os.path import join

path_to_summary_table = join(os.environ['SPECTRASTACKS_DIR'], "results", "table_v0.data")

lineWavelength,Survey,Redshift,L_MIN,L_MAX,L_MEAN,N_in_stack,R_stack,spm_light_age,spm_light_age_err_plus,spm_light_age_err_minus,spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus,spm_stellar_mass,spm_stellar_mass_err_plus,spm_stellar_mass_err_minus,spm_EBV,gp_EBV_4862_4341,gp_EBV_4862_4341_err,gp_EBV_4862_4102,gp_EBV_4862_4102_err,gp_BD_4102_4341,gp_BD_4102_4341_err,gp_SFR_O2_3728,gp_SFR_O2_3728_err,gp_SFR_H1_4862,gp_SFR_H1_4862_err,gp_12logOH_tremonti04,gp_12logOH_tremonti04_err,gp_12logOH_tremonti04_intrinsic,gp_12logOH_tremonti04_intrinsic_err = n.loadtxt(path_to_summary_table)	


# EBV comparison
spm_EBV
gp_EBV_4862_4341,gp_EBV_4862_4341_err,

spm_EBV
gp_EBV_4862_4102,gp_EBV_4862_4102_err

# metallicity comparison
spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus
gp_12logOH_tremonti04,gp_12logOH_tremonti04_err,

spm_light_metallicity,spm_light_metallicity_err_plus,spm_light_metallicity_err_minus
gp_12logOH_tremonti04_intrinsic,gp_12logOH_tremonti04_intrinsic_err

# mass metallicity relation

x_obs_1 = spm_stellar_mass - 0.32 * gp_SFR_O2_3728 
x_obs_2 = spm_stellar_mass - 0.32 * gp_SFR_H1_4862 
y_obs_1 = spm_light_metallicity
y_obs_2 = gp_12logOH_tremonti04
y_obs_3 = gp_12logOH_tremonti04_intrinsic

x_rel = n.arange(8.5,11.5,0.05) - 10 
y_pr = y_rel(x_rel)
y_rel = lambda x_rel : 8.90 + 0.39*x_rel − 0.20*x_rel*x_rel − 0.077*x_rel*x_rel*x_rel + 0.064*x_rel*x_rel*x_rel*x_rel


p.ylabel(r'$12 \log(O/H)$')
p.xlabel('$\log(M_*)-0.32 \log(SFR)-10$')
p.xlim((9,11.5))
p.savefig()
