import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import numpy as n
import os

import sys

import ClusterScalingRelations as clsr
from scipy.interpolate import interp1d
import StellarMass as sm
smhmr = sm.StellarMass()
scl = clsr.ClusterScalingRelations_Mantz2016()

cat = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'catCluster-SPIDERS_RASS_CLUS-v2.0.fits'))[1].data
spm = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster_statistics_2016-11-08-DR14_spm.fits'))[1].data

#previous catalogs
m2x = []
lx = []
mass = []
zs = []
for cc in cat:
	gal = (spm['CLUS_ID']==cc['CLUS_ID'])&(spm['Chabrier_ELODIE_stellar_mass']>0.)
	#&(abs(spm['Z']-cc['SCREEN_CLUZSPEC'])<0.2)
	mhs = n.arange(7,16,0.01)
	itp = interp1d( scl.logM500_to_L( mhs ,cc['SCREEN_CLUZSPEC']), mhs)
	itMS = interp1d( smhmr.SMHMr(10**mhs,cc['SCREEN_CLUZSPEC'])*10**mhs, mhs)
	mass_to_X = itMS(spm['Chabrier_ELODIE_stellar_mass'][gal]) -itp(cc['LX0124'])
	ooo = n.ones_like(spm['Chabrier_ELODIE_stellar_mass'][gal])
	m2x.append(mass_to_X)
	lx.append(ooo*cc['LX0124'])
	mass.append(spm['Chabrier_ELODIE_stellar_mass'][gal])
	zs.append(ooo*cc['SCREEN_CLUZSPEC'])
	
	


logm2x = n.hstack((m2x))	
ok = (n.isnan(logm2x)==False)&(logm2x != -n.inf)&(logm2x != n.inf)
bins=n.arange(-7, 0.5, 0.1)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])

x = n.hstack((lx))
y = n.hstack((mass))
z = n.hstack((zs))

p.figure(1, (5,5))
p.title('SPIDERS DR14 galaxies')
p.scatter(x, y, s=8, c=z, edgecolor='none')
cb = p.colorbar()
cb.set_label('redshift')
p.xlabel('LX [erg/s]')
p.ylabel('stellar mass [Msun]')
p.xscale('log')
p.yscale('log')
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'LX-mass.png'))
p.clf()



p.figure(1, (5,5))
p.title('SPIDERS DR14 galaxies')
p.plot(spm['Z'], spm['Chabrier_ELODIE_stellar_mass'], 'b	,', alpha=0.5, label='targets')
p.plot(z, y, 'r,', label='cluster members')
p.xlabel('redshift')
p.ylabel('stellar mass [Msun]')
#p.xscale('log')
p.yscale('log')
p.xlim((0,0.7))
p.ylim((1e9,1e12))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'redshift-mass.png'))
p.clf()
