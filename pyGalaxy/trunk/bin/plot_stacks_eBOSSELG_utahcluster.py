import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile

from scipy.stats import chi2

import sys
import os 
from os.path import join
import glob
eb17_sL = n.array(glob.glob("elg270_eboss17*Z1*_stac*.fits"))
eb67_sL = n.array(glob.glob("elg270_eboss67*Z1*_stac*.fits"))
eb67_sL.sort()
eb17_sL.sort()

for ii in range(len(eb17_sL)):
	hd_eb17 = fits.open(eb17_sL[ii])
	fig = p.figure(1,(11,6))
	fig.add_subplot(211)
	p.plot(hd_eb17[1].data['wavelength'],hd_eb17[1].data['meanStack'],'k',label=eb17_sL[ii].split('_')[-3]+"<z<"+eb17_sL[ii].split('_')[-2])
	p.ylabel(r'$f_\lambda \; 10^{-17}$ erg/cm2/s/A')
	p.xlim((2250,7000))
	p.ylim((0.05,2.))
	p.yscale('log')
	p.legend(loc=2)
	p.grid()
	p.title('eboss17')
	fig.add_subplot(212)
	p.plot(hd_eb17[1].data['wavelength'],hd_eb17[1].data['meanStack']/hd_eb17[1].data['jackknifStackErrors'],'k')
	p.xlabel('wavelength [A]')
	p.ylabel('SNR per pixel')
	p.xlim((2250,7000))
	p.ylim((1,1000))
	p.yscale('log')
	p.legend(loc=2)
	p.savefig(eb17_sL[ii][:-5]+'.png')
	p.clf()


for ii in range(len(eb67_sL)):
	hd_eb67 = fits.open(eb67_sL[ii])
	fig = p.figure(1,(11,6))
	fig.add_subplot(211)
	p.plot(hd_eb67[1].data['wavelength'],hd_eb67[1].data['meanStack'],'k',label=eb67_sL[ii].split('_')[-3]+"<z<"+eb67_sL[ii].split('_')[-2])
	p.ylabel(r'$f_\lambda \; 10^{-17}$ erg/cm2/s/A')
	p.xlim((2250,7000))
	p.ylim((0.05,2.))
	p.yscale('log')
	p.legend(loc=2)
	p.grid()
	p.title('eboss67')
	fig.add_subplot(212)
	p.plot(hd_eb67[1].data['wavelength'],hd_eb67[1].data['meanStack']/hd_eb67[1].data['jackknifStackErrors'],'k')
	p.xlabel('wavelength [A]')
	p.ylabel('SNR per pixel')
	p.xlim((2250,7000))
	p.ylim((1,1000))
	p.yscale('log')
	p.legend(loc=2)
	p.savefig(eb67_sL[ii][:-5]+'.png')
	p.clf()



hdu_eb17 = fits.open("elg270_eboss17summaryTable_stack_comparison.fits")
hdu_eb67 = fits.open("elg270_eboss67summaryTable_stack_comparison.fits")


fig = p.figure(1,(11,4.5))
fig.add_subplot(121)
ccp=p.scatter(hdu_eb17[1].data['Z_1'], hdu_eb17[1].data['chi2_Z1'], s=5, c=hdu_eb17[1].data['gmag'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
p.colorbar(shrink=0.7)
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss17')

fig.add_subplot(122)
ccp=p.scatter(hdu_eb67[1].data['Z_1'], hdu_eb67[1].data['chi2_Z1'], s=5, c=hdu_eb67[1].data['gmag'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
cb=p.colorbar(shrink=0.7)
cb.set_label('g')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss67')
p.legend()

p.savefig('chi2-z-elg-v591-gMag.png')
p.clf()


fig = p.figure(1,(11,4.5))
fig.add_subplot(121)
ccp=p.scatter(hdu_eb17[1].data['Z_1'], hdu_eb17[1].data['chi2_Z1'], s=5, c=hdu_eb17[1].data['grcol'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
p.colorbar(shrink=0.7)
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss17')

fig.add_subplot(122)
ccp=p.scatter(hdu_eb67[1].data['Z_1'], hdu_eb67[1].data['chi2_Z1'], s=5, c=hdu_eb67[1].data['grcol'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
cb=p.colorbar(shrink=0.7)
cb.set_label('g-r')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss67')
p.legend()

p.savefig('chi2-z-elg-v591-grColor.png')
p.clf()

fig = p.figure(1,(11,4.5))
fig.add_subplot(121)
ccp=p.scatter(hdu_eb17[1].data['Z_1'], hdu_eb17[1].data['chi2_Z1'], s=5, c=hdu_eb17[1].data['rzcol'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
p.colorbar(shrink=0.7)
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss17')

fig.add_subplot(122)
ccp=p.scatter(hdu_eb67[1].data['Z_1'], hdu_eb67[1].data['chi2_Z1'], s=5, c=hdu_eb67[1].data['rzcol'], marker = 'o', rasterized= True)
ccp.set_edgecolors('face')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
cb=p.colorbar(shrink=0.7)
cb.set_label('r-z')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.title('eboss67')
p.legend()

p.savefig('chi2-z-elg-v591-rzColor.png')
p.clf()




fig = p.figure(0,(9,4.5))
fig.add_subplot(121)
p.plot(hdu_eb17[1].data['Z_1'], hdu_eb17[1].data['chi2_Z1'], 'k+', label='eboss17', rasterized= True)
x=n.arange(0,2.5,0.01)
y=0.75+(x-0.85)**2/3.
p.plot(x,y, 'b--')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

fig.add_subplot(122)
p.plot(hdu_eb67[1].data['Z_1'], hdu_eb67[1].data['chi2_Z1'], 'k+', label='eboss67', rasterized= True)
p.plot(x,y, 'b--')
p.xlabel('best redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

p.savefig('chi2-z-elg-v591.png')
p.clf()

sel = (hdu_eb17[1].data['Z_1']<0.5) | (hdu_eb17[1].data['Z_1']>1.2)
n.median(hdu_eb17[1].data['chi2_Z2'][sel] - hdu_eb17[1].data['chi2_Z1'][sel])
n.median(hdu_eb17[1].data['chi2_Z3'][sel] - hdu_eb17[1].data['chi2_Z1'][sel])

fig = p.figure(0,(9,4.5))
fig.add_subplot(121)
p.plot(hdu_eb17[1].data['Z_2'][sel], hdu_eb17[1].data['chi2_Z2'][sel], 'k+', label='eboss17 z1<0.5 or z1>1.2', rasterized= True)
p.plot(x,y, 'b--')
p.xlabel('second redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

sel = (hdu_eb67[1].data['Z_1']<0.5) | (hdu_eb67[1].data['Z_1']>1.2)

fig.add_subplot(122)
p.plot(hdu_eb67[1].data['Z_2'][sel], hdu_eb67[1].data['chi2_Z2'][sel], 'k+', label='eboss67 z1<0.5 or z1>1.2', rasterized= True)
p.plot(x,y, 'b--')
p.xlabel('second redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

p.savefig('chi2-z-elg-v591_Z2.png')
p.clf()

sel = (hdu_eb17[1].data['Z_1']<0.5) | (hdu_eb17[1].data['Z_1']>1.2)

fig = p.figure(0,(9,4.5))
fig.add_subplot(121)
p.plot(hdu_eb17[1].data['Z_3'][sel], hdu_eb17[1].data['chi2_Z3'][sel], 'k+', label='eboss17 z1<0.5 or z1>1.2', rasterized= True)
p.plot(x,y, 'b--')
p.xlabel('third redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

sel = (hdu_eb67[1].data['Z_1']<0.5) | (hdu_eb67[1].data['Z_1']>1.2)

fig.add_subplot(122)
p.plot(hdu_eb67[1].data['Z_3'][sel], hdu_eb67[1].data['chi2_Z3'][sel], 'k+', label='eboss67 z1<0.5 or z1>1.2', rasterized= True)
p.plot(x,y, 'b--')
p.xlabel('third redshift')
p.ylabel(r'$\chi^2/dof$')
p.xlim((0,2.5))
p.ylim((0.7,1.5))
p.grid()
p.legend()

p.savefig('chi2-z-elg-v591_Z3.png')
p.clf()
