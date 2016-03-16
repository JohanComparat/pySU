import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize

#saundersFct=lambda v, A, logv0,alpha,sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( -n.log10( 1 + v /10**logv0)/(2*sig**2.))
#schechterFct=lambda v, A, logv0,alpha, sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( - v / 10**logv0 /(2*sig**2.) )
#ps * (10**logl/10**logls)**(a+1) * n.e**(-10**logl/10**logls)
#doublePL=lambda v,A,logv0 ,a,b: 10**A * 10**( (1+a)*( v - 10**logv0) + b )
#vf0 = lambda v, A, v0, alpha, beta : n.log10( 10**A * (1+ (10**v/10**v0)**beta) * n.e**(- (10**v/10**v0)**alpha ) )
vf = lambda v, A, v0, alpha, beta, s0 : n.log10( 10**A * (1+ (10**v/10**v0)**beta) * n.e**( -n.log10( 1 +10**v/10**(v0 ))**alpha / (2*(s0 )**alpha)))

# limits at z0
p0 = n.array([-4.3, 2.9, 2.6, -2.8, 0.17])
p0S = n.array([-3.3, 2.5, 1.7, -2.6, 0.12])
error = 0.05

limits_04 = [70, 3000]
limits_10 = [200, 3000]
limits_25 = [400, 3000]
limits_40 = [650, 5000]
zmin = 0.
zmax = 5.

NDecimal = 3

# defining directories :
dir = join("D:","data","MultiDark")
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

dir_boxes =  n.array([dir_04, dir_10, dir_25, dir_40])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])
qty_limits = n.array([limits_04, limits_10, limits_25, limits_40])
volume_boxes =  n.array([400.**3., 1000**3., 2500**3., 4000.**3.])

property_dir = "vmax-mvir"
type = "hist"
cos = "Central" #"Central" # centrak or satellite ?
qty = "vmax"


f=open(join("outputs_vmax", "outputs-vmax-cen-fits.txt"),'w')

f.write( "we consider the "+type+" of "+qty+" of "+ cos+ "\n")
f.write( "in the redshift range"+str(zmin)+" "+str(zmax)+ "\n")
#f.write( "for the boxes"+dir_boxes+ "\n")
#f.write( zList_files
f.write( "within the following limits for each box"+str(qty_limits) +"\n")
f.write( "each box has a volume of"+str(volume_boxes)+ "Mpc3/h3" +"\n")


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Satellite-vmax_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Satellite-vmax_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Satellite-vmax_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Satellite-vmax_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= 0.0) & (z_04 <= zmax)
s_10 = (z_10 >= 0.0) & (z_10 <= zmax)
s_25 = (z_25 >= 0.0) & (z_25 <= zmax)
s_40 = (z_40 >= 0.0) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
f.write("all redshifts available:"+str(set(redshift)) + "\n")
vmax = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
f.write( "min and max masses available:"+str( n.min(vmax))+" "+str( n.max(vmax))+ "\n")
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
f.write("min and max Y available:"+str(n.min(yData))+" "+ str(n.max(yData))+ "\n")
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
f.write("min and max Y error available:"+str( n.min(yDataErr))+ " "+str( n.max(yDataErr)) + "\n")


################################

def fitSingleRedshiftM200Function(z0,z1):
	print z0, z1
	s_04 = (z_04 < z1) & (z_04 > z0 )
	s_10 = (z_10 < z1) & (z_10 > z0 )
	s_25 = (z_25 < z1) & (z_25 > z0 )
	s_40 = (z_40 < z1) & (z_40 > z0 )

	redshift = n.hstack(( z_04[s_04], z_10[s_10] ))#, z_25[s_25], z_40[s_40]))
	f.write( "all redshifts available:"+str(set(redshift))+"\n")
	vmax = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10])) )#, xData_25[s_25], xData_40[s_40])))
	f.write( "min and max masses available:"+str(n.min(vmax))+" "+str(n.max(vmax))+"\n")
	yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10])))#, yData_25[s_25], yData_40[s_40])))
	f.write( "min and max Y available:"+str( n.min(yData))+" "+str( n.max(yData))+"\n")
	yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10])))#, yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
	f.write( "min and max Y error available:"+str( n.min(yDataErr))+" "+str( n.max(yDataErr))+"\n")

	# with curve fit
	f.write( "with curve fit \n")
	popt, cov = curve_fit(vf, vmax, yData,sigma=error*yData, p0 = p0 , maxfev = 5000000)
	A0, vcut0, a0, b0, s0 = n.round(popt,NDecimal)

	f.write( "redshift 0 model for the SAT vmax cumulative function :"+"\n")
	f.write( "A(z=0) & = "+str(A0)+ r"\pm "+str( n.round(cov[0][0]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" M_{vir}^{cut}(z=0) & = "+str(vcut0)+ r"\pm "+str( n.round(cov[1][1]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" \alpha(z=0) & = "+str(a0)+ r"\pm "+str( n.round(cov[2][2]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" \beta(z=0) & = "+str(b0)+ r"\pm "+str( n.round(cov[3][3]**0.5, NDecimal))+ r'\\'+"\n")

	gg=open("outputs_vmax\MF-sat-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
	n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
	gg.close()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
	p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
	#p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
	#p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
	xModel = n.arange(n.min(vmax),n.max(vmax),0.1)
	yModel_CF = vf(xModel, A0, vcut0, a0, b0,s0)
	p.plot(xModel, yModel_CF,'k--',label="model")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n$_{sat}$(>M)') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0)+" "+str(s0))
	p.grid()
	p.savefig("outputs_vmax\MF-sat-4params-z-"+str(n.mean(redshift))+"-fit.png")
	p.clf()
	
	yModel_CF = vf(vmax, A0, vcut0, a0, b0,s0)

	p.figure(1,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(vmax,10**yData/10**yModel_CF,'bo')
	p.axhline(1)
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n$_{sat}(>M)$ data / model') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.ylim((.95,1.05))
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
	p.grid()
	p.savefig("outputs_vmax\MF-sat-4params-z-"+str(n.mean(redshift))+"dataovermodel.png")
	p.clf()


fitSingleRedshiftM200Function(-0.1,0.01)
fitSingleRedshiftM200Function(0.091,0.15)
fitSingleRedshiftM200Function(0.25,0.35)
fitSingleRedshiftM200Function(0.65,0.7)
fitSingleRedshiftM200Function(0.99,1.1)
fitSingleRedshiftM200Function(1.4,1.5)
fitSingleRedshiftM200Function(2.,2.2)
fitSingleRedshiftM200Function(2.8,3.)


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Central-vmax_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Central-vmax_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Central-vmax_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Central-vmax_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= 0.0) & (z_04 <= zmax)
s_10 = (z_10 >= 0.0) & (z_10 <= zmax)
s_25 = (z_25 >= 0.0) & (z_25 <= zmax)
s_40 = (z_40 >= 0.0) & (z_40 <= zmax)


################################ Plot cumulative halo mass function and model at z=0  ################################

def fitSingleRedshiftM200Function(z0,z1):

	s_04 = (z_04 < z1) & (z_04 > z0 )
	s_10 = (z_10 < z1) & (z_10 > z0 )
	s_25 = (z_25 < z1) & (z_25 > z0 )
	s_40 = (z_40 < z1) & (z_40 > z0 )

	redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
	f.write( "all redshifts available:"+str( set(redshift))+"\n")
	vmax = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
	f.write( "min and max masses available:"+str( n.min(vmax))+" "+str( n.max(vmax))+"\n")
	yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
	f.write( "min and max Y available:"+str( n.min(yData))+" "+str( n.max(yData))+"\n")
	yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
	f.write( "min and max Y error available:"+str(n.min(yDataErr))+" "+str( n.max(yDataErr))+"\n")

	# with curve fit
	f.write( "with curve fit"+"\n")
	popt, cov = curve_fit(vf, vmax, yData,sigma=error*yData, p0 = p0 , maxfev = 5000000)
	A0, vcut0, a0, b0,s0 = n.round(popt,NDecimal)

	f.write( "redshift 0 model for the vmax cumulative function :"+"\n")
	f.write( "A(z=0) & = "+str(A0)+ r"\pm "+str( n.round(cov[0][0]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" M_{vir}^{cut}(z=0) & = "+str(vcut0)+ r"\pm "+str( n.round(cov[1][1]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" \alpha(z=0) & = "+str(a0)+ r"\pm "+str( n.round(cov[2][2]**0.5, NDecimal))+ r'\\'+"\n")
	f.write( r" \beta(z=0) & = "+str(b0)+ r"\pm "+str( n.round(cov[3][3]**0.5, NDecimal))+ r'\\'+"\n")

	gg=open("outputs_vmax\MF-cen-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
	n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
	gg.close()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
	p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
	p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
	p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
	xModel = n.arange(n.min(vmax),n.max(vmax),0.1)
	yModel_CF = vf(xModel, A0, vcut0, a0, b0,s0)
	p.plot(xModel, yModel_CF,'k--',label="model")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n$_{cen}$(>M)') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0)+" "+str(s0))
	p.savefig("outputs_vmax\MF-cen-4params-z-"+str(n.mean(redshift))+"-fit.png")
	p.clf()
	
	yModel_CF = vf(vmax, A0, vcut0, a0, b0, s0)

	p.figure(1,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(vmax,10**yData/10**yModel_CF,'bo')
	p.axhline(1)

	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n$_{cen}$(>M) data / model') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.ylim((.95,1.05))
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
	p.grid()
	p.savefig("outputs_vmax\MF-cen-4params-z-"+str(n.mean(redshift))+"dataovermodel.png")
	p.clf()


fitSingleRedshiftM200Function(-0.1,0.01)
fitSingleRedshiftM200Function(0.091,0.15)
fitSingleRedshiftM200Function(0.25,0.35)
fitSingleRedshiftM200Function(0.65,0.7)
fitSingleRedshiftM200Function(0.99,1.1)
fitSingleRedshiftM200Function(1.4,1.5)
fitSingleRedshiftM200Function(2.,2.2)
fitSingleRedshiftM200Function(2.8,3.)


f.close()

sys.exit()

