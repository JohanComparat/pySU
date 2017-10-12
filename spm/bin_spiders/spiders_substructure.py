import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits
import astropy.units as u

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

cat = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'validatedclusters_catalogue_2016-07-04-DR14_version_round1-v4_Xmass-v1.fits.gz'))[1].data
spm = fits.open(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'validatedclusters_catalogue_2016-07-04-DR14_version_round1-v4_Xmass-v1_spm.fits'))[1].data

volume_rough = aa.comoving_volume(0.5)*2200.*n.pi/129600
volume = volume_rough.value

#previous catalogs
m2x = []
lx = []
mass = []
zs = []
m200c = []
crd = []
for cc in cat:
	gal = (spm['CLUS_ID']==cc['CLUS_ID'])
	#&(spm['Chabrier_MILES_stellar_mass']>0.)
	member = (cc['SCREEN_ISMEMBER_W'][:len(gal.nonzero()[0])]==1)
	#&(abs(spm['Z']-cc['SCREEN_CLUZSPEC'])<0.2)
	
	stellar_mass = spm['Chabrier_MILES_stellar_mass'][gal][member]
	sm =stellar_mass[stellar_mass>0]
	
	mhs = n.arange(7,16,0.01)
	itp = interp1d( scl.logM500_to_L( mhs ,cc['SCREEN_CLUZSPEC']), mhs)
	
	itMS = interp1d( smhmr.SMHMr(10**mhs,cc['SCREEN_CLUZSPEC'])*10**mhs, mhs)
	
	mass_to_X = itMS(sm) - n.log10(cc['M200c']) #itp(cc['LX0124'])
	ooo = n.ones_like(sm)
	m2x.append(mass_to_X)
	lx.append(ooo*cc['LX0124'])
	m200c.append(ooo*cc['M200c'])
	mass.append(sm)
	zs.append(ooo*cc['SCREEN_CLUZSPEC'])
	crd.append(aa.critical_density(cc['SCREEN_CLUZSPEC']).to(u.solMass/u.megaparsec**3.).value)

x = n.hstack((lx))
y = n.hstack((mass))
z = n.hstack((zs))
m200c = n.hstack((m200c))
crd = n.hstack((crd))

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
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'LX-mass.png'))
p.clf()

p.figure(1, (5,5))
p.title('SPIDERS DR14 galaxies')
p.plot(spm['Z'], spm["Chabrier_MILES_stellar_mass"], 'b,', label='targets')
p.plot(z, y, 'r,', label='cluster members')
p.xlabel('redshift')
p.ylabel('stellar mass [Msun]')
#p.xscale('log')
p.yscale('log')
p.xlim((0,0.7))
p.ylim((1e9,1e12))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster',  'redshift-mass.png'))
p.clf()

logm2x = n.hstack((m2x))	
bins=n.arange(-7, 0.5, 0.1)

basis = (n.isnan(logm2x)==False)&(logm2x != -n.inf)&(logm2x != n.inf)

arbitrary_factor =5.

p.figure(1, (5,5))
ok = (basis)&(x>1e44)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>44')

ok = (basis)&(x>10**44.5)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>44.5')

ok = (basis)&(x>1e45)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='LX>45')

ok = (basis)&(m200c>10**14)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='M200c>14', ls='dashed')

ok = (basis)&(m200c>10**15)
out = n.log10(n.histogram(logm2x[ok], bins=bins)[0])
p.plot((bins[1:]+bins[:-1])/2., n.log10(out/arbitrary_factor), label='M200c>15', ls='dashed')

xs = n.arange(-7, 0.01, 0.01)
logfsat= lambda logxi, a, b, logN0, exponent : n.log10( 10**logN0 * (10**logxi)**a)# * n.e**(-b*(10**logxi)**exponent))

p.plot(xs, logfsat(xs, -0.81, 5.81, -2.25, -2.54), label='-0.81')
p.plot(xs, logfsat(xs, -0.18, 5.81, -1.2, -.54), label='-0.18')

p.xlabel('log10(SMHMR(stellar mass) / HaloMass(Lx ray))')
p.ylabel('histogram')
#p.xscale('log')
#p.yscale('log')
p.ylim((-1.5, 0.5))
p.xlim((-4,0))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'LX-mass-histogram.png'))
p.clf()

