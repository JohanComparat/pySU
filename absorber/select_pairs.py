"""
To do a first estimate the number of pairs (with ra, dec ), you can use the following two files :
- clusters: https://www.sdss.org/dr16/data_access/value-added-catalogs/?vac_id=spiders-x-ray-galaxy-cluster-catalogue-for-dr16 where ra, dec, redshift are given by RA, DEC, SCREEN_CLUZSPEC. R200C_DEG is the Apparent R200c radius of the galaxy cluster, I would recommend to record matches up to 2 or 3 x R200c_DEG. (Indeed the cluster center is not well constrained: RA, DEC are not very precise)
- all SDSS spectra: https://firefly.mpe.mpg.de/v1_1_0/v5_13_0/catalogs/spAll-v5_13_0.fits  where ra, dec, redshift are given by PLUG_RA, PLUG_DEC, Z
i.e. something like :
& distance (PLUG_RA, PLUG_DEC // RA, DEC ) < 3 x R200C_DEG
& SCREEN_CLUZSPEC + 0.01 < Z

Then among the pairs, filter the ones with the correct redshift combination: cluster -- background object i.e. for each possible line obtain the total number of pairs with the right redshift separation (SCREEN_CLUZSPEC , Z )

Depending on how many pairs of cluster-QSO are in the data, one could try to look at trends as a function of cluster-centric distance using either X-ray center (RA, DEC) of optical center (RA_OPT, DEC_OPT).

"""
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from matplotlib.patches import Polygon
from astropy import constants
import astropy.io.fits as fits
from matplotlib import gridspec
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from sklearn.neighbors import BallTree
from astropy.table import Table, Column

import pymangle 
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
#
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
match_radius_arcmin = interp1d( np.arange(0.001,1.5,0.001), cosmo.arcsec_per_kpc_comoving(np.arange(0.001,1.5,0.001))*1200*u.kpc/(60*u.arcsec) )

deg_to_rad = np.pi / 180.
nl = lambda selection : len(selection.nonzero()[0])

# define pathes
env = 'HOME'
p2_spall = os.path.join(os.environ[env], 'data2/firefly/v1_1_0/v5_13_0/catalogs/spAll-v5_13_0.fits')
p_2_out = os.path.join(os.environ[env], 'wwwDir/stuff/catalogue_qso_v5_13_0.fits')

p2_spall = os.path.join(os.environ[env], 'data2/firefly/v1_1_0/26/catalogs/specObj-SDSS-dr12.fits')
p_2_out = os.path.join(os.environ[env], 'wwwDir/stuff/catalogue_qso_26.fits')

p2_codex_bcg = os.path.join(os.environ[env], 'hegcl/SPIDERS/mastercatalogue_FINAL_CODEXID.fits')


spall = fits.open(p2_spall)#[1].data
codex = fits.open(p2_codex_bcg)[1].data

#ok = (spall[1].data['Z']>0.3) & (spall[1].data['Z_ERR'] > 0 ) & (spall[1].data['Z_ERR']<spall[1].data['Z']) & (spall[1].data['ZWARNING']==0) # & (spall[1].data['CLASS']=="QSO")
ok = (spall[1].data['Z']>0.3) & (spall[1].data['Z_ERR'] > 0 ) & (spall[1].data['Z_ERR']<spall[1].data['Z']) & (spall[1].data['ZWARNING']==0) & (spall[1].data['CLASS']=="QSO")
dr16_rough = ( spall[1].data['PLUG_DEC'] < 15 ) & (spall[1].data['PLUG_RA']>100 ) & (spall[1].data['PLUG_RA']<280 )
select = (dr16_rough==False) & (ok)
plug_ra  = spall[1].data['PLUG_RA'][select]
plug_dec = spall[1].data['PLUG_DEC'][select]
Z_dr16 = spall[1].data['Z'][select]
print(len(Z_dr16))

plate   = spall[1].data['PLATE'][select]
mjd     = spall[1].data['MJD'][select]
fiberid = spall[1].data['FIBERID'][select]


#clu_coord = deg_to_rad * np.array([codex['DEC_OPT'], codex['RA_OPT']]).T
#Tree_Cluster_Cat = BallTree(clu_coord, metric='haversine')

CAT_coord = deg_to_rad * np.transpose([plug_dec, plug_ra])
CAT_Tree = BallTree(CAT_coord, metric='haversine')

# around each cluster, query 
def get_data(id_c = 0, N_RADIUS = 3):
	clu_coord = deg_to_rad * np.transpose([np.array([codex['DEC_OPT'][id_c]]), np.array([codex['RA_OPT'][id_c]])])
	radius = deg_to_rad * codex['R200C_DEG'][id_c]

	indexes, distances = CAT_Tree.query_radius(clu_coord, r = radius * N_RADIUS, return_distance = True, sort_results = True)
	relevant = (Z_dr16[indexes[0]] > codex['SCREEN_CLUZSPEC'][id_c] + 0.01)

	Z_out        = Z_dr16[indexes[0][relevant]]
	distance_out = distances[0][relevant] / deg_to_rad
	plate_out    = plate   [indexes[0][relevant]]
	mjd_out      = mjd     [indexes[0][relevant]]
	fiberid_out  = fiberid [indexes[0][relevant]]
	plug_ra_out  = plug_ra [indexes[0][relevant]]
	plug_dec_out = plug_dec[indexes[0][relevant]]

	cluster_ID          = np.ones_like(Z_out).astype('str')
	cluster_ID[:] = codex['CLUS_ID'][id_c]
	cluster_ra          = codex['RA_OPT'][id_c]  * np.ones_like(Z_out) 
	cluster_dec         = codex['DEC_OPT'][id_c] * np.ones_like(Z_out)
	cluster_z           = codex['SCREEN_CLUZSPEC'][id_c] * np.ones_like(Z_out)
	cluster_r200c_deg   = codex['R200C_DEG'][id_c] * np.ones_like(Z_out)
	cluster_KT          = codex['KT'][id_c] * np.ones_like(Z_out)

	DATA_i = np.transpose([ 
		cluster_ID       , 
		cluster_ra       , 
		cluster_dec      , 
		cluster_z        , 
		cluster_r200c_deg, 
		cluster_KT       , 
		Z_out            ,
		distance_out     ,
		plate_out        ,
		mjd_out          ,
		fiberid_out      ,
		plug_ra_out      ,
		plug_dec_out     
		])

	print(id_c, DATA_i.shape)
	return DATA_i

d0 = get_data(id_c = 0)
for id_c in np.arange(1, len(codex), 1) :
	d1 = get_data(id_c)
	d0 = np.vstack((d0,d1))

t = Table()

t.add_column( Column(name="cluster_ID"        , data = d0.T[0] ) )
t.add_column( Column(name="cluster_ra"        , unit='deg', data = d0.T[1].astype('float') ) )
t.add_column( Column(name="cluster_dec"       , unit='deg', data = d0.T[2].astype('float') ) )
t.add_column( Column(name="cluster_z"         , data = d0.T[3].astype('float') ) )
t.add_column( Column(name="cluster_r200c_deg" , data = d0.T[4].astype('float') ) )
t.add_column( Column(name="cluster_KT"        , data = d0.T[5].astype('float') ) )
t.add_column( Column(name="galaxy_z"          , data = d0.T[6].astype('float') ) )
t.add_column( Column(name="angular_separation", unit='deg'      , data = d0.T[7].astype('float') ) )
t.add_column( Column(name="plate"         , data = d0.T[8].astype('int') ) )
t.add_column( Column(name="mjd"           , data = d0.T[9].astype('int') ) )
t.add_column( Column(name="fiberid"       , data = d0.T[10].astype('int') ) )
t.add_column( Column(name="RA"       , data = d0.T[11].astype('float') ) )
t.add_column( Column(name="DEC"       , data = d0.T[12].astype('float') ) )

t.write(p_2_out, overwrite=True)

