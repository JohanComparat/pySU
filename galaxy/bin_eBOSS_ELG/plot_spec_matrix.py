import numpy as n
import glob
import os, sys
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n

from cycler import cycler
# 1. Setting prop cycle on default rc parameter
p.rc('lines', linewidth=1.3)
p.rc('axes', prop_cycle=cycler('color', ["#638bd9", "#e586b6", "#dd7c25"]) )#, "#598664", "#9631a9", "#cff159", "#534f55", "#ab1519", "#89dbde"]) )

path_2_fly_dir = os.path.join(os.environ['HOME'], 'wwwDir/sdss/elg/stacks/v2/firefly/')
path_2_spec_dir = os.path.join(os.environ['HOME'], 'SDSS/stacks/v4/')

file_list = n.array(glob.glob(os.path.join(path_2_spec_dir,"*O2EW*specMatrix.dat")))
file_list.sort()
baseNames = n.array([ os.path.basename(bn) for bn in file_list ])

ii=0

path_2_MAT = os.path.join( path_2_spec_dir, baseNames[ii] )
path_2_spec = os.path.join( path_2_spec_dir, baseNames[ii][:-15] )
path_2_asc = os.path.join( path_2_spec_dir, baseNames[ii][:-21]+'.asc' )
path_2_figure = os.path.join( path_2_fly_dir, 'specMatrix-'+baseNames[ii]+'.png')
title_2 = baseNames[ii][10:-21]
title_2_sp = title_2.split('_')
title = str(n.round(float(title_2_sp[0]),2))+'<z<'+str(n.round(float(title_2_sp[2]),2))+', '+str(n.round(float(title_2_sp[3]),2))+'<'+title_2_sp[4]+'<'+str(n.round(float(title_2_sp[5]),2))


#def plot_it(path_2_spec, path_2_spFly, path_2_figure, title):
full_list = n.loadtxt(path_2_asc, unpack=True)
redshift = full_list[-1]

print(path_2_figure)
matrix = n.loadtxt(path_2_MAT)
path_2_spec)

d = fits.open(path_2_spec)
sel = (d[1].data['NspectraPerPixel']>0.5*n.max(d[1].data['NspectraPerPixel']))&(d[1].data['medianStack']>0)
x_data = d[1].data['wavelength'][sel]
y_data = d[1].data['medianStack'][sel]
y_err = d[1].data['jackknifStackErrors'][sel]
N_spec = d[1].data['NspectraPerPixel'][sel]

wavelength = d[1].data['wavelength']

dLambda = 0.0001
l_start=3.25 
l_end=3.6
wave= 10**n.arange(l_start, l_end, dLambda)


# A4 figure
fig = p.figure(0, (10.2, 12.2), frameon=False )

# panel top left: spectrum + models
fig.add_subplot(4,1,1, 
				title=title, 
				#xlabel='wavelength [Angstrom, rest frame]', 
				ylabel=r'Flux [$f_\lambda 10^{-17}$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]',
				ylim=((n.min(y_data)*0.9,1.1* scoreatpercentile(y_data,90))),
				xlim=((XMIN, XMAX)) )

p.plot(x_data, y_data, color='grey', lw=3, label='stack, N ELG='+str(int(n.min(N_spec))) )
p.plot(x_model, y_model, label='firefly MILE UV extended model')
p.legend(loc=0, fontsize=12, frameon=False)
p.grid()

## panel top right: residuals
fig.add_subplot(4,1,2, 
			ylabel='|(data-model)|/error', 
			xlabel='wavelength [Angstrom, rest frame]', 
			xlim=((XMIN, XMAX)), 
			yscale='log')#, xlim=((-2.5,2.5)))

s_model = (x_model > XMIN) & (x_model < XMAX)
s_data  = (x_data  > XMIN) & (x_data  < XMAX)

itp_model = interp1d( x_model[s_model], y_model[s_model] )
itp_data  = interp1d( x_data[s_data], y_data[s_data] )

chis = abs(y_data[s_data]-itp_model(x_data[s_data]))/y_err[s_data]

p.plot(x_data[s_data], chis)
##p.legend( frameon=False)#, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
p.grid()

### age metal
### log10 light weighted age [Gyr] 
fig.add_subplot(4,4,9, xlabel=r'Z [$\log_{10}(Z/Z_\odot)$]', ylabel=r'age  [$\log_{10}(age/yr)$]')
x  = n.array([ m[2].header['metallicity_massW']    ])
xe = n.array([ m[2].header['metallicity_massW_up_1sig'], m[2].header['metallicity_massW_low_1sig']       ])
y  = n.log10(10**9 * 10**n.array([ m[2].header['age_massW']    ])  )
ye = n.log10(10**9 * 10**n.array([ m[2].header['age_massW_up_1sig'], m[2].header['age_massW_low_1sig']       ]) )
p.plot(n.array([x, x]), ye, lw=1, color='b')

x  = n.array([ m[2].header['metallicity_massW']    ])
xe = n.array([ m[2].header['metallicity_massW_up_1sig'], m[2].header['metallicity_massW_low_1sig']       ])
y  = n.log10(10**9 * 10**n.array([ m[2].header['age_massW']    ])  )
ye = n.log10(10**9 * 10**n.array([ m[2].header['age_massW_up_1sig'], m[2].header['age_massW_low_1sig']       ]) )
p.plot(x,y,'bo')
p.plot(xe, n.array([y, y]), lw=1, color='b')

p.grid()

## mass, age
fig.add_subplot(4,4,10, xlabel=r'mass [$\log_{10}(M/M_\odot)$]') # , ylabel='age  [$\log_{10}(age/yr)$]'
y  = n.log10(10**9 * 10**n.array([ m[2].header['age_massW']    ])  )
ye = n.log10(10**9 * 10**n.array([ m[2].header['age_massW_up_1sig'], m[2].header['age_massW_low_1sig']       ]) )
x  = n.array([ m[2].header['total_mass'] ])
xe = n.array([ m[2].header['total_mass_up_1sig'], m[2].header['total_mass_low_1sig'] ])
p.plot(n.array([x, x]), ye, lw=1, color='b')

y  = n.log10(10**9 * 10**n.array([ m[2].header['age_massW']    ])  )
ye = n.log10(10**9 * 10**n.array([ m[2].header['age_massW_up_1sig'], m[2].header['age_massW_low_1sig']       ]) )
x  = n.array([ m[2].header['total_mass'] ])
xe = n.array([ m[2].header['total_mass_up_1sig'], m[2].header['total_mass_low_1sig'] ])
p.plot(x,y,'bo')
p.plot(xe, n.array([y, y]), lw=1, color='b')

p.grid()

## mass metal
fig.add_subplot(4,4,13, xlabel=r'Z [$\log_{10}(Z/Z_\odot)$]', ylabel=r'mass [$\log_{10}(M/M_\odot)$]')
y  = n.array([ m[2].header['total_mass'] ])
ye = n.array([ m[2].header['total_mass_up_1sig'], m[2].header['total_mass_low_1sig'] ])
x  = n.array([ m[2].header['metallicity_massW']    ])
xe = n.array([ m[2].header['metallicity_massW_up_1sig'], m[2].header['metallicity_massW_low_1sig']       ])
p.plot(n.array([x, x]), ye, lw=1, color='b')

y  = n.array([ m[2].header['total_mass'] ])
ye = n.array([ m[2].header['total_mass_up_1sig'], m[2].header['total_mass_low_1sig'] ])
x  = n.array([ m[2].header['metallicity_massW']    ])
xe = n.array([ m[2].header['metallicity_massW_up_1sig'], m[2].header['metallicity_massW_low_1sig']       ])
p.plot(x,y,'bo')
p.plot(xe, n.array([y, y]), lw=1, color='b')

p.grid()


## mass EBV
fig.add_subplot(4,4,14, xlabel=r'E(B-V)')#, ylabel='mass [$\log_{10}(M/M_\odot)$]')
y  = n.array([ m[2].header['total_mass'] ])
x  = n.array([ m[2].header['EBV']    ])
ye = n.array([ m[2].header['total_mass_up_1sig'], m[2].header['total_mass_low_1sig'] ])
p.plot(n.array([x, x]), ye, lw=1)

p.grid()
## weigths age
fig.add_subplot(4,2,6, xlabel=r'age  [$\log_{10}(age/yr)$]', ylabel=r'mass weight SSP', ylim=(-0.02,1.02) )
weightM = n.array([ m[2].header['weightMass_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
##weightL = n.array([ m[2].header['weightLight_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
age = n.array([ n.log10(10**9 * 10**m[2].header['log_age_ssp_'+str(jj)]) for jj in n.arange(m[2].header['ssp_number'])])
##metal = n.array([ m[2].header['metal_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
iiid=n.argsort(age)
p.plot(age[iiid], weightM[iiid], marker='+', ls='', markersize=12, color='b')

p.grid()

## weigths metallicity
fig.add_subplot(4,2,8	, xlabel=r'Z [$\log_{10}(Z/Z_\odot)$]', ylabel=r'mass weight SSP', ylim=(-0.02,1.02) )
weightM = n.array([ m[2].header['weightMass_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
##weightL = n.array([ m[2].header['weightLight_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
##age = n.array([ n.log10(10**9 * 10**m[2].header['log_age_ssp_'+str(jj)]) for jj in n.arange(m[2].header['ssp_number'])])
metal = n.array([ m[2].header['metal_ssp_'+str(jj)] for jj in n.arange(m[2].header['ssp_number'])])
iiid=n.argsort(metal)
p.plot(metal[iiid], weightM[iiid], marker='+', ls='', markersize=12, color='b')

p.grid()
p.tight_layout()
p.savefig(path_2_figure)
p.clf()
