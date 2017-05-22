import sys
import os
import numpy as n
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 2000000
matplotlib.rcParams.update({'font.size': 13})
#matplotlib.use('Agg')
import matplotlib.pyplot as p
from cycler import cycler
from scipy.spatial import KDTree
from scipy.stats import norm

# ADD FILTER ON CHI2/NDOF

"""
Plots stellar mass error and stellar mass vs. mass and redshift
"""

z_bins = n.arange(0, 4.1, 0.1)
m_bins = n.arange(8.5,12.6,1.)
sn_bins = n.array([0, 0.5, 1, 2, 10, 100]) # n.logspace(-1.,2,6)

err_bins = n.arange(-2., 2.2, 0.1) 
x_err_bins = (err_bins[1:] + err_bins[:-1])/2.

#n.hstack(( n.array([-10000., -10., -5.]), n.arange(-2.5, 2.5, 0.1), n.array([5., 10., 10000.]) )) 

#prefix = 'SDSS'
#hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxySdss26.fits'))
#redshift_reliable = (hdus[1].data['Z'] >= 0) & ( hdus[1].data['Z_ERR'] >= 0) & (hdus[1].data['ZWARNING'] == 0) & (hdus[1].data['Z'] > hdus[1].data['Z_ERR'] )

prefix = 'BOSS'
out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results', 'catalogs')

hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxyEbossDR14.fits'))


NN, bb = n.histogram(hdus[1].data['PLATE'], bins = n.arange(0, 11000,1))
repeat_plates = bb[:-1][(NN>4000)]

rp = repeat_plates[4]
sel_rp = (hdus[1].data['PLATE'] == rp)
data = hdus[1].data[sel_rp]
tree = KDTree(n.transpose([hdus[1].data['PLUG_RA'][sel_rp], hdus[1].data['PLUG_DEC'][sel_rp]]))
ids = n.array(tree.query_ball_tree(tree, 0.0001))
mjds = n.array(list(set(hdus[1].data['MJD'][sel_rp])))
mjd_sel = (hdus[1].data['MJD'][sel_rp] == mjds[-1])

id_sort = ids[mjd_sel]

imf = 'Chabrier'
diff = []
for jj in range(len(id_sort)):
	ids = list(n.copy(id_sort[jj]))
	id_highest_snr = n.argmax(data[ids]['SN_MEDIAN_ALL'])
	del ids[id_highest_snr]

	M_1   = 10**data[imf+'_stellar_mass'][id_highest_snr]
	M_1_e = M_1 * (10**abs((data[imf+'_stellar_mass' + '_err_plus'][id_highest_snr] - data[imf+'_stellar_mass' + '_err_minus'][id_highest_snr])/2.)-1.)

	M_2s   = 10**data[imf+'_stellar_mass'][ids]
	M_2_es = M_1 * (10**abs((data[imf+'_stellar_mass' + '_err_plus'][ids] - data[imf+'_stellar_mass' + '_err_minus'][ids])/2.)-1.)
	z_agree = (data['ZWARNING_NOQSO'][id_highest_snr]==0)&(abs(data['Z_NOQSO'][id_highest_snr] - data['Z_NOQSO'][ids]) < 0.05)
	if M_1 > M_1_e and M_1 > 0. and data['Z_ERR_NOQSO'][id_highest_snr]>0 and data['Z_NOQSO'][id_highest_snr]> data['Z_ERR_NOQSO'][id_highest_snr] : 
		arr = (M_1 - M_2s)/(M_1_e**2. + M_2_es**2.)**0.5
		#print arr[z_agree]
		diff.append( arr[z_agree] )

ps = [-0.15, -0.2, 0.05, 0.25, 0.65, 0.35]

x_norm = n.arange(-2,2,0.01) 

out, xxx = n.histogram(n.hstack((diff)), bins = err_bins, normed=True)
outN = n.histogram(n.hstack((diff)), bins = err_bins)[0]
N100 = (outN>10)

eb = p.errorbar(x_err_bins[N100], out[N100], xerr=(xxx[1:][N100]-xxx[:-1][N100])/2., yerr = out[N100]*outN[N100]**(-0.5), fmt='+', color='r')

p.plot(x_norm, ps[4]*norm.pdf(x_norm, loc=ps[0], scale=ps[2]), label='N('+str(ps[0])+','+str(ps[2])+')', ls='dashed', lw=0.5)
p.plot(x_norm, ps[5]*norm.pdf(x_norm, loc=ps[1], scale=ps[3]), label='N('+str(ps[1])+','+str(ps[3])+')', ls='dashed', lw=0.5)
p.plot(x_norm, ps[5]*norm.pdf(x_norm, loc=ps[1], scale=ps[3]) + ps[4]*norm.pdf(x_norm, loc=ps[0], scale=ps[2]), lw=0.5)
p.yscale('log')
p.ylabel('pdf')
p.xlim((-1.1, 1.1))
p.ylim((1e-2, 10.))
p.grid()
p.xlabel(r'$(M_1-M_2)/\sqrt{\sigma^2_{M1}+\sigma^2_{M2}}$')
p.savefig(os.path.join(out_dir, "pdf_diff_mass_repeat.jpg" ))
p.show()
