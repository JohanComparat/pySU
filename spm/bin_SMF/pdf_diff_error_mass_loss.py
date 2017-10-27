from lib_spm import *

from cycler import cycler

from scipy.stats import norm


# ADD FILTER ON CHI2/NDOF

"""
Plots stellar mass error and stellar mass vs. mass and redshift
"""

z_bins = n.arange(0, 4.1, 0.1)
m_bins = n.arange(7.5,12.6,2.)
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

redshift_reliable =  (hdus[1].data['CLASS_NOQSO'] == "GALAXY") & (hdus[1].data['Z_NOQSO'] >= 0) & ( hdus[1].data['Z_ERR_NOQSO'] >= 0) & (hdus[1].data['ZWARNING_NOQSO'] == 0) & (hdus[1].data['Z_NOQSO'] > hdus[1].data['Z_ERR_NOQSO'] ) # (hdus[1].data['SN_MEDIAN_ALL'] > 0.1 ) & 

error_reliable = (hdus[1].data['Chabrier_stellar_mass_err_plus'] > hdus[1].data['Chabrier_stellar_mass_err_minus'] ) & (hdus[1].data['Salpeter_stellar_mass_err_plus'] > hdus[1].data['Salpeter_stellar_mass_err_minus'] )  & (hdus[1].data['Kroupa_stellar_mass_err_plus'] > hdus[1].data['Kroupa_stellar_mass_err_minus'] ) & (hdus[1].data['Chabrier_stellar_mass_err_plus'] > 0. ) & ( hdus[1].data['Chabrier_stellar_mass_err_minus'] > 0. ) & (hdus[1].data['Salpeter_stellar_mass_err_plus'] > 0. ) & ( hdus[1].data['Salpeter_stellar_mass_err_minus'] > 0. )  & (hdus[1].data['Kroupa_stellar_mass_err_plus'] > 0. ) & ( hdus[1].data['Kroupa_stellar_mass_err_minus'] > 0. ) & (hdus[1].data['Chabrier_stellar_mass_err_plus'] < 10. ) & ( hdus[1].data['Chabrier_stellar_mass_err_minus'] < 10. ) & (hdus[1].data['Salpeter_stellar_mass_err_plus'] < 10. ) & ( hdus[1].data['Salpeter_stellar_mass_err_minus'] < 10. )  & (hdus[1].data['Kroupa_stellar_mass_err_plus'] < 10. ) & ( hdus[1].data['Kroupa_stellar_mass_err_minus'] < 10. )

mass_reliable = (hdus[1].data['Chabrier_stellar_mass'] > 0 ) & ( hdus[1].data['Chabrier_stellar_mass'] < 14. ) & (hdus[1].data['Salpeter_stellar_mass'] > 0 ) & ( hdus[1].data['Salpeter_stellar_mass'] < 14. )  & (hdus[1].data['Kroupa_stellar_mass'] > 0 ) & ( hdus[1].data['Kroupa_stellar_mass'] < 14. ) 

ok = (error_reliable) & (mass_reliable) & (redshift_reliable)

snr = hdus[1].data['SN_MEDIAN_ALL'][ok]
zz = hdus[1].data['Z_NOQSO'][ok]

M_c   = 10**hdus[1].data['Chabrier_stellar_mass'][ok]
M_c_e = M_c * (10**abs((hdus[1].data['Chabrier_stellar_mass' + '_err_plus'][ok] - hdus[1].data['Chabrier_stellar_mass' + '_err_minus'][ok])/2.) -1)

M_s   = 10**hdus[1].data['Salpeter_stellar_mass'][ok]
M_s_e = M_s * (10**abs((hdus[1].data['Salpeter_stellar_mass' + '_err_plus'][ok] - hdus[1].data['Salpeter_stellar_mass' + '_err_minus'][ok])/2.) -1)

M_k   = 10**hdus[1].data['Kroupa_stellar_mass'][ok]
M_k_e = M_k * (10**abs((hdus[1].data['Kroupa_stellar_mass' + '_err_plus'][ok] - hdus[1].data['Kroupa_stellar_mass' + '_err_minus'][ok])/2.) -1)


def plot_errPDF(mass, dMs, title, ps):
	p.figure(1, (4.5, 4.5))
	p.axes([0.2,0.2,0.7,0.7])
	cc=cycler('color', ['r', 'g', 'b', 'm', 'k'])

	for ii, (mb, c) in enumerate(zip(m_bins[:-1],cc)):
		out, xxx = n.histogram(dMS[(mass>m_bins[ii])&(mass<m_bins[ii+1])], bins = err_bins, normed=True)
		outN = n.histogram(dMS[(mass>m_bins[ii])&(mass<m_bins[ii+1])], bins = err_bins)[0]
		N100 = (outN>10)
		if len(x_err_bins[N100])>10 :
			eb = p.errorbar(x_err_bins[N100], out[N100], xerr=(xxx[1:][N100]-xxx[:-1][N100])/2., yerr = out[N100]*outN[N100]**(-0.5), label=str(mb)+"-"+str(m_bins[ii+1]), fmt='+', color=c['color'])

	x_norm = n.arange(-2,2,0.01) 
	p.plot(x_norm, norm.pdf(x_norm, loc=0, scale=1), label='N(0,1)', ls='solid', lw=0.5)
	p.plot(x_norm, ps[4]*norm.pdf(x_norm, loc=ps[0], scale=ps[2]), label='N('+str(ps[0])+','+str(ps[2])+')', ls='dashed', lw=0.5)
	p.plot(x_norm, ps[5]*norm.pdf(x_norm, loc=ps[1], scale=ps[3]), label='N('+str(ps[1])+','+str(ps[3])+')', ls='dashed', lw=0.5)
	p.plot(x_norm, ps[5]*norm.pdf(x_norm, loc=ps[1], scale=ps[3]) + ps[4]*norm.pdf(x_norm, loc=ps[0], scale=ps[2]), lw=0.5)
	p.legend(loc=2, frameon=False)
	p.yscale('log')
	p.ylabel('pdf')
	p.xlim((-1.5, 1.5))
	p.ylim((1e-2, 10.))
	p.grid()
	p.title(title)
	p.xlabel(r'$(M_1-M_2)/\sqrt{\sigma^2_{M1}+\sigma^2_{M2}}$')
	p.savefig(os.path.join(out_dir, "pdf_diff_mass_"+title+".jpg" ))
	p.clf()


diff = (M_c - M_s)/(M_c_e**2. + M_s_e**2.)**0.5
dMS=diff
ps = [-0.2, -0.25, 0.1, 0.5, 0.5, 0.5]
plot_errPDF(n.log10(M_c), dMS, title='Chabrier-Salpeter', ps=ps)

diff = (M_c - M_k)/(M_c_e**2. + M_k_e**2.)**0.5
dMS=diff
ps = [0., 0.05, 0.05, 0.4, 0.7, 0.3]
plot_errPDF(n.log10(M_c), dMS, title='Chabrier-Kroupa',ps=ps)

diff = (M_s - M_k)/(M_s_e**2. + M_k_e**2.)**0.5
dMS=diff
ps = [+0.15, +0.3, 0.1, 0.4, 0.5, 0.5]
plot_errPDF(n.log10(M_s), dMS, title='Salpeter-Kroupa',ps=ps)


