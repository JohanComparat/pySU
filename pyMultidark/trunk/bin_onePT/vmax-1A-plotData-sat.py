import glob
import sys
import cPickle
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.interpolate import interp1d

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "vmax", "MD_vmax_summary.fits") )[1].data
boxes = set(data['boxName'])
mk = {"MD_0.4Gpc": '1', "MD_1Gpc": 'x' ,"MD_2.5Gpc": '2',"MD_4Gpc": '3' ,"MD_2.5GpcNW": '4', "MD_4GpcNW": '+'}

NminCount = 1 #1000
nolim = [0,4000]
limits_04 = nolim # [100, 400]
limits_10 = nolim # [250, 1000]
limits_25 = nolim # [600, 1300]
limits_40 = nolim # [1200, 1600]
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.001

def plotData(qty = 'vmax', cos = "sat", zmin = -0.01, zmax = 2.3):
	"""
	Plots the data to be used in the fits later in the analysis.
	"""
	# redshift selection
	zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
	# mass selection
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>limits_25[0]) &(data["log_"+qty+"_max"]<limits_25[1])) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>limits_40[0])&(data["log_"+qty+"_max"]<limits_40[1]))
	# minimum number counts selection
	nSel = (data['dN_counts_'+cos]>NminCount)
	# altogether
	ok = (zSel) & (mSel) & (nSel)
	
	# x coordinates
	log_vmax = (n.log10(data["log_"+qty+"_min"][ok])+n.log10(data["log_"+qty+"_max"][ok]))/2.
	vmax = 10**log_vmax
	print len(vmax), n.min(vmax), n.max(vmax)
	
	# y coordinates
	norm = (100)**3. /(cosmo.H(data["redshift"][ok]).value)**6.
	log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos][ok])
	
	log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"][ok])
	print n.min(log_VF), n.max(log_VF)
	
	MD04=(data["boxName"]=="MD_0.4Gpc")
	MD10=(data["boxName"]=="MD_1Gpc")
	MD25=(data["boxName"]=="MD_2.5Gpc")
	MD25N=(data["boxName"]=="MD_2.5GpcNW")
	MD40=(data["boxName"]=="MD_4Gpc")
	MD40N=(data["boxName"]=="MD_4GpcNW")
	
	# now the plots
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"][ok]==box)
		#print box, mk[box]
		sc1=p.scatter(log_vmax[bb], log_VF[bb], c=data["redshift"][ok][bb], s=20, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')
	
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'log$_{10} [(V^3/H^3(z)\; dn(V)/dV]$') # log$_{10}[ 
	gl = p.legend(loc=3,fontsize=12)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	# p.xlim((1.5, 3.5))
	#p.ylim((-3.5,-1))
	p.grid()
	p.title('satellite data')
	p.savefig(join(dir,qty,"vmax-"+cos+"-differential-function-data.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"][ok]==box)
		#print box, mk[box]
		sc1=p.scatter(log_vmax[bb], log_VF_c[bb], c=data["redshift"][ok][bb], s=30, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')

	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'log$_{10} [(V^3/H^3(z)\; n(>V)]$') # log$_{10}[ 
	gl = p.legend(loc=3,fontsize=12)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	# p.xlim((1.5, 3.5))
	#p.yscale('log')
	p.grid()
	p.title('satellite data')
	p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-data.png"))
	p.clf()
	"""
	log_vmax_04 = (n.log10(data["log_"+qty+"_min"][ok & MD04])+n.log10(data["log_"+qty+"_max"][ok & MD04]))/2.
	error_04 = data['std90_pc_cen'][ok & MD04]
	
	log_vmax_10 = (n.log10(data["log_"+qty+"_min"][ok & MD10])+n.log10(data["log_"+qty+"_max"][ok & MD10]))/2.
	error_10 = data['std90_pc_cen'][ok & MD10]

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_vmax_04, 100*error_04, c=data["redshift"][ok & MD04], s=5, marker='o',label="MD04", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	# p.xlim((1.5, 3.5))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"vmax-"+cos+"-data04-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_vmax_10, 100*error_10, c=data["redshift"][ok & MD10], s=5, marker='o',label="MD10", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	# p.xlim((1.5, 3.5))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"vmax-"+cos+"-data10-uncertainty.png"))
	p.clf()
	"""

plotData(qty = 'vmax', cos = 'sat')
