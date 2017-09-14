import os
import numpy as n
import idlsave
import glob
import astropy.io.fits as fits

spm_dir = os.path.join(os.environ['DATA_DIR'], "spm")
gama_dir = os.path.join(spm_dir, "GAMAmock")
input_dir = os.path.join(gama_dir, "inputs")

RecycleFraction = 0.43
h0 = 0.693
mass_factor = 1. #/(h0*(1.-RecycleFraction))

cat = fits.open(os.path.join(gama_dir, "catalogs","GAMA.mock.spectra.spm.fits"))
out_name = os.path.join(gama_dir, "catalogs","GAMA.mock.spectra.spm.in.out.fits")

file_list = n.array(glob.glob(os.path.join(input_dir, 'mock_input_GAMA', 'input_*_GAMA_M10_z0.15.idl')))
#file_list.sort()

result_list =n.array(glob.glob(os.path.join(gama_dir, 'results', 'gal_*_GAMA_M10_z0.15.fits')))
#result_list.sort()

def get_in_values(identifier, dV=-9999.99):
	file_i = os.path.join(input_dir, 'mock_input_GAMA', 'input_'+identifier+'_GAMA_M10_z0.15.idl')
	#identifier = os.path.basename(file_i).split('_')[1]
	file_r = os.path.join(gama_dir, 'results', 'gal_'+identifier+'_GAMA_M10_z0.15.fits')
	print( os.path.basename(file_i), os.path.basename(file_r) )
	if os.path.isfile(file_i) and os.path.isfile(file_r):
		spec = idlsave.read(file_i)
		model = fits.open(file_r)
		# gets integrated quantities from the input
		stellar_mass = spec['sfh_bulge'].sum()/mass_factor + spec['sfh_disk'].sum()/mass_factor 
		age_massW = n.log10(n.sum(10**spec['ages_lg']* (spec['sfh_bulge'] + spec['sfh_disk']))/stellar_mass)
		
		z_massW_bulge = n.sum(spec['z_bulge'][(spec['sfh_bulge']>0)]* spec['sfh_bulge'][(spec['sfh_bulge']>0)])/n.sum(spec['sfh_bulge'][(spec['sfh_bulge']>0)]) 
		z_massW_disk = n.sum(spec['z_disk'][(spec['sfh_disk']>0)] * spec['sfh_disk'][(spec['sfh_disk']>0)])/n.sum(spec['sfh_disk'][(spec['sfh_disk']>0)])
		if z_massW_disk>0 and z_massW_bulge>0 :
			z_massW = z_massW_disk + z_massW_bulge
		elif z_massW_disk>0 and n.isnan(z_massW_bulge):
			z_massW = z_massW_disk
		elif z_massW_bulge>0 and n.isnan(z_massW_disk):
			z_massW = z_massW_bulge
			stellar_mass
		return stellar_mass, age_massW, z_massW
	else :
		return dV, dV, dV

stellar_mass = n.zeros(len(cat[1].data))
age_massW = n.zeros(len(cat[1].data))
z_massW = n.zeros(len(cat[1].data))

for jj, index in enumerate(cat[1].data['ID']):
	identifier = str(int(index)).zfill(4)
	stellar_mass[jj], age_massW[jj], z_massW[jj] = get_in_values(identifier)

# now write the merged catalog
c1 = fits.Column(name='input_stellar_mass', format='D', unit='log10(M_sun)', array=n.log10(stellar_mass))
c2 = fits.Column(name='input_age_massW', format='D', unit='Gyr', array=10**age_massW)
c3 = fits.Column(name='input_metallicity_massW', format='D', unit='Z_sun', array=z_massW)


new_columns = cat[1].data.columns + c1 + c2 + c3 
hdu = fits.BinTableHDU.from_columns(new_columns)

if os.path.isfile(out_name):
	os.remove(out_name)

hdu.writeto(out_name)

#coldefs = fits.ColDefs([c1, c2])
#tbhdu = fits.BinTableHDU.from_columns(coldefs)

#

# compares input and corresponding output
#print spec['ages_lg']
#print model[1].header['age_ssp_0'], model[1].header['age_ssp_1']
#
#print model[1].header['stellar_mass']
#print model[1].header['age_massW'], age_massW
#