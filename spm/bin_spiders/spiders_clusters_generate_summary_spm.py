import astropy.io.fits as fits
import os
import time
import numpy as n
t0 = time.time()

dir = os.path.join(os.environ['DATA_DIR'], 'spiders', 'cluster')

hd = fits.open(os.path.join(dir, 'validatedclusters_catalogue_2016-07-04-DR14_version_round1-v4_Xmass-v1.fits.gz'))[1].data

spm = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'sdss_eboss_firefly-dr14.fits'))
spm_plate = spm[1].data['PLATE']
spm_mjd = spm[1].data['MJD']
spm_fiberid = spm[1].data['FIBERID']

index = n.arange(spm[1].header['NAXIS2'])

ids = []
id_clus = []
for el in hd :
	print( el['CLUS_ID'] )
	for plate, mjd, fib in zip(	el['ALLPLATE'], el['ALLMJD'], el['ALLFIBERID'] ):
		if plate != -99999 :
			print( plate, mjd, fib )
			correspondence = (spm_plate==plate )&( spm_mjd == mjd )&( spm_fiberid == fib )
			iid = index[correspondence]
			print( iid, time.time()-t0 )
			ids.append(iid)
			id_clus.append(el['CLUS_ID'])

id_save = ids
ids = n.hstack((ids))
DATA = spm[1].data[ids]
import astropy.io.fits as fits
import os
import time
import numpy as n
t0 = time.time()

all_cols = []
for cc in DATA.columns:
	all_cols.append(fits.Column(name = cc.name, format = cc.format, array=DATA[cc.name]))
	print( len(DATA[cc.name]))
# add a CLUS_ID column

col0 = fits.Column(name='CLUS_ID',format='A10', array=n.array(id_clus) )
all_cols.append(col0)

#define the table hdu 
hdu_cols  = fits.ColDefs(all_cols)
tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
#define the header
prihdr = fits.Header()
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])
thdulist.writeto(os.path.join(dir, 'cluster_statistics_2016-11-08-DR14_spm.fits'))
	