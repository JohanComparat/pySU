import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import angles
#import AngularSeparation
from astropy import coordinates as coord

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

# get cluster center

# distance to center
# rescale to r200c_deg
# get the latest min(ages) of the ssp
# compute SFR

# now looks at individual galaxies
# and gets the highest SFR for each galaxy
# youngest age
highest_sfrs = []
youngest_ages = []
sep_r200c = []

for cc in cat:
	center = coord.ICRS(ra=cc['RA_OPT']*u.degree, dec=cc['DEC_OPT']*u.degree)
	gal = (spm['CLUS_ID']==cc['CLUS_ID'])
	#all_members = coord.ICRS()
	#separations = center.separation(all_members)/(cc['R200C_DEG']*u.degree)).value
	for id_cc, (pla, mjd, fib) in enumerate(zip(cc['ALLPLATE'][:len(gal.nonzero()[0])], cc['ALLMJD'][:len(gal.nonzero()[0])], cc['ALLFIBERID'][:len(gal.nonzero()[0])])):
		sel = (gal) & (spm['PLATE']==pla) & (spm['MJD']==mjd) & (spm['FIBERID']==fib)
		if  len(sel.nonzero()[0])>0 :
			n_cp = spm['Chabrier_MILES_nComponentsSSP'][sel].astype('int')[0]
			if n_cp > 0 :
				all_ages = n.array([ spm['Chabrier_MILES_age_ssp_'+str(ii)][sel][0] for ii in n.arange(n_cp) ])
				all_masses = n.array([ spm['Chabrier_MILES_stellar_mass_ssp_'+str(ii)][sel][0] for ii in n.arange(n_cp) ])
				sfr_inst = all_masses / all_ages
				youngest_ages.append(n.min(all_ages))
				highest_sfrs.append(n.max(sfr_inst))
				position = coord.ICRS(cc['ALLRA'][id_cc]*u.degree, cc['ALLDEC'][id_cc]*u.degree)
				sep_r200c.append( (center.separation(position)/(cc['R200C_DEG']*u.degree)).value )


highest_sfrs = n.array(highest_sfrs)
youngest_ages = n.array(youngest_ages)
sep_r200c = n.array(sep_r200c)

p.figure(1, (5,5))
p.title('SPIDERS')
p.plot(sep_r200c, highest_sfrs, 'r+')
p.xlabel('r/r200c')
p.ylabel('SFR [Msun/yr]')
#p.xscale('log')
p.yscale('log')
p.xlim((0.08,1.5))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'disteance-2-center-SFR.png'))
p.clf()

dx = ( n.max(sep_r200c) - n.min(sep_r200c) ) /3.
r_b = n.arange(n.min(sep_r200c), n.max(sep_r200c) + dx, dx)

p.figure(1, (5,5))
for ii,bb in enumerate(r_b[:-1]):
	sub = (sep_r200c>bb)&(sep_r200c<r_b[ii+1])
	p.hist(highest_sfrs[sub], label=str(n.round(bb,3))+"<"+str(n.round(r_b[ii+1],3)), cumulative=True, normed=True, histtype='step')

p.ylabel('normed cumulative distribution')
p.xlabel('SFR [Msun/yr]')
p.xscale('log')
p.ylim((-0.01, 1.01))
p.grid()
p.legend(frameon=False, loc=0)
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'disteance-2-center-SFR-histograms.png'))
p.clf()


p.figure(1, (5,5))
p.title('SPIDERS')
p.plot(sep_r200c, youngest_ages, 'r+')
p.xlabel('r/r200c')
p.ylabel('age [yr]')
p.xscale('log')
p.yscale('log')
p.xlim((0.1,5))
p.grid()
p.savefig(os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster', 'disteance-2-center-AGE.png'))
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

