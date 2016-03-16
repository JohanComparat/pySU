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

# limits at z0
limits_04 = [1e10, 5e12]
limits_10 = [5e11, 5e13]
limits_25 = [5e12, 5e14]
limits_40 = [1e13, 5e15]
zmin = 0.
zmax = 5

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

property_dir = "M200c-mvir"
type = "hist"
cos = "Central" #"Central" # centrak or satellite ?
qty = "M200c"

f=open(join("outputs", "outputs-m200c-sat-fits-all-at-once.txt"),'w')

f.write( "we consider the "+type+" of "+qty+" of "+ cos+ "\n")
f.write( "in the redshift range"+str(zmin)+" "+str(zmax)+ "\n")
#f.write( "for the boxes"+dir_boxes+ "\n")
#f.write( zList_files
f.write( "within the following limits for each box"+str(qty_limits) +"\n")
f.write( "each box has a volume of"+str(volume_boxes)+ "Mpc3/h3" +"\n")

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= 0.0) & (z_04 <= zmax)
s_10 = (z_10 >= 0.0) & (z_10 <= zmax)
s_25 = (z_25 >= 0.0) & (z_25 <= zmax)
s_40 = (z_40 >= 0.0) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10] ))#, z_25[s_25], z_40[s_40]))
f.write( "all redshifts available:"+str(set(redshift))+"\n")
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10])) )#, xData_25[s_25], xData_40[s_40])))
f.write( "min and max masses available:"+str(n.min(M200c))+" "+str(n.max(M200c))+"\n")
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10])))#, yData_25[s_25], yData_40[s_40])))
f.write( "min and max Y available:"+str( n.min(yData))+" "+str( n.max(yData))+"\n")
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10])))#, yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
f.write( "min and max Y error available:"+str( n.min(yDataErr))+" "+str( n.max(yDataErr))+"\n")

#first guess

vcut0 = 12.826 # 13.8425
vcut1 = 1.0 # -0.8762
vcut2 =0.0 #  0.0414

b0 = -0.9 # -0.8749
b1 = 0.0 # -0.0124
b2 = 0.0 # -0.0090

a0 = 0.578# 0.5702
a1 = -0.2 # -0.072
a2 = 0.0 # 0.014

A0 = -3.72 # -4.0178
A1 = 0.9 # 0.8232
A2 = 0.0 # -0.0714


vfG = lambda v, z, A0, A1, vcut0, vcut1, a0, a1, b0, b1 : n.log10( 10**(A0 + A1 * z) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 +a1*z) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7])
chi2fun = lambda ps : n.sum( (vfGbis(M200c,redshift,ps) - yData)**2. / (yData*0.017)**2. )/(len(yData) - 8)
p1 = n.array([ A0, A1, vcut0, vcut1, a0, a1, b0, b1 ])

f.write( "looks for the optimum parameters"+"\n")
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

pOpt = res.x
cov = res.direc
chi2perpoint = lambda ps : (vfGbis(M200c,redshift,ps) - yData)**2. / (yData*0.017)**2. 
chi2pp = chi2perpoint(pOpt)

f.write( "A(z) & = "+str(n.round(pOpt[0],NDecimal))+ " + "+str(n.round(pOpt[1],NDecimal))+r'\times z \\'+"\n")
f.write( r" M_{cut}(z) & = "+str(n.round(pOpt[2],NDecimal))+"+ "+str(n.round(pOpt[3],NDecimal))+r'\times z \\'+"\n")
f.write( r" \alpha(z) & = "+str(n.round(pOpt[4],NDecimal))+" + "+str(n.round(pOpt[5],NDecimal))+r'\times z \\'+"\n") #+ '+str(a2)+r'\times z^2 \\'
f.write( r" \beta(z) & = "+str(n.round(pOpt[6],NDecimal))+" + "+str(n.round(pOpt[7],NDecimal))+r'\times z \\'+"\n")


#f.write( "A(z) & = ("+str(n.round(pOpt[0],NDecimal))+ r"\pm "+str( n.round(cov[0][0]**0.5, NDecimal))+") + ("+str(n.round(pOpt[1],NDecimal))+ r"\pm "+str( n.round(cov[1][1]**0.5, NDecimal))+r')\times z \\'+"\n")
#f.write( r" M_{cut}(z) & = ("+str(n.round(pOpt[2],NDecimal))+ r"\pm "+str( n.round(cov[2][2]**0.5, NDecimal))+") + ("+str(n.round(pOpt[3],NDecimal))+ r"\pm "+str( n.round(cov[3][3]**0.5, NDecimal))+r')\times z \\'+"\n")
#f.write( r" \alpha(z) & = ("+str(n.round(pOpt[4],NDecimal))+ r"\pm "+str( n.round(cov[4][4]**0.5, NDecimal))+") + ("+str(n.round(pOpt[5],NDecimal))+ r"\pm "+str( n.round(cov[5][5]**0.5, NDecimal))+r')\times z \\'+"\n") #+ '+str(a2)+r'\times z^2 \\'
#f.write( r" \beta(z) & = ("+str(n.round(pOpt[6],NDecimal))+ r"\pm "+str( n.round(cov[6][6]**0.5, NDecimal))+") + ("+str(n.round(pOpt[7],NDecimal))+ r"\pm "+str( n.round(cov[7][7]**0.5, NDecimal))+r')\times z \\'+"\n")

xModel = n.arange(n.min(M200c),14,0.1)
X,Y = n.meshgrid(xModel,n.arange(0,n.max(redshift)+0.025,0.025))
Z = vfGbis(X,Y,pOpt)

n.savetxt(join("outputs","M200c-sat-cumulative-function-best_fit.txt"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )

#######################################################
# now plots the results of the fit
print "now plots the results of the fit"

vmax_mod, z_mod, n_mod = n.loadtxt(join("outputs","M200c-sat-cumulative-function-best_fit.txt"), unpack=True)

tpl=(n_mod>-6)
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(vmax_mod[tpl], n_mod[tpl], c=z_mod[tpl],s=5, marker='o',label="model", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.ylim((-7,0))
p.xlim((9.5, 14))
p.grid()
p.savefig(join("outputs","M200c-sat-cumulative-function-model.png"))
p.clf()


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(M200c, yData, c=redshift,s=5, marker='o',label="data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.ylim((-7,0))
p.xlim((9.5, 14))
p.grid()
p.savefig(join("outputs","M200c-sat-cumulative-function-data.png"))
p.clf()



p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(M200c, redshift, c=chi2pp,s=5, marker='o', rasterized=True, vmin = 0, vmax = 4)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("$\chi^2$")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r'z') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.xlim((9.5, 14))
p.grid()
p.savefig(join("outputs","M200c-sat-cumulative-function-chi2pp.png"))
p.clf()



sys.exit()























sys.exit()

vfG = lambda v, z, A1, vcut1, a1, b1 : n.log10( 10**(A0 + A1 * z) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 +a1*z) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1],ps[2],ps[3])
chi2fun = lambda ps : n.sum( (vfGbis(M200c,redshift,ps) - yData)**2. / (n.ones_like(yData)*0.05)**2. )/len(yData)
p1 = n.array([ 0., 0., 0., 0.])
print "looks for the optimum parameters"
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

print "ndof=",len(yData)
print res
A1, vcut1, a1, b1 = n.round(res.x,4)

print "A(z) & = "+str(A0)+" + "+str(A1)+r'\times z \\'
print r" M_{cut}(z) & = "+str(vcut0)+" + "+str(vcut1)+r'\times z \\'
print r" \alpha(z) & = "+str(a0)+" + "+str(a1)+r'\times z \\' #+ '+str(a2)+r'\times z^2 \\'
print r" \beta(z) & = "+str(b0)+" + "+str(b1)+r'\times z \\'


vfG = lambda v, z, vcut1, b1 : n.log10( 10**(A0 ) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 ) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1])
chi2fun = lambda ps : n.sum( (vfGbis(M200c,redshift,ps) - yData)**2. / (n.ones_like(yData)*0.05)**2. )/len(yData)

p1 = n.array([ 0.,  0.])

print "looks for the optimum parameters"
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

print "ndof=",len(yData)
print res
A1, vcut1, a1, b1 = n.round(res.x,4)

print "A(z) & = "+str(A0)+" + "+str(A1)+r'\times z \\'
print r" M_{cut}(z) & = "+str(vcut0)+" + "+str(vcut1)+r'\times z \\'
print r" \alpha(z) & = "+str(a0)+" + "+str(a1)+r'\times z \\' #+ '+str(a2)+r'\times z^2 \\'
print r" \beta(z) & = "+str(b0)+" + "+str(b1)+r'\times z \\'


#####################

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
dd = ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)
sc1 = ax.scatter(M200c, redshift,yData, s=n.ones_like(yData)*3, c='r', marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
ax.legend()
ax.set_xlabel(r'log $M_{200c}$ [km s$^{-1}$]')
ax.set_ylabel('redshift')
ax.set_ylim((0,n.max(redshift)))
ax.set_zlabel(r'log N($>M_{200c}$) [ h$^3$ Mpc$^{-3}$]')
ax.set_zlim((-8,0))
p.show()

#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir , "M200c-cumulative-function-allZ-model.pdf"))
p.clf()

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
sc1 = ax.scatter(M200c, redshift, 10**yData/10**vfGbis(M200c,redshift,res.x), s=n.ones_like(yData)*3, c='r', marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
ax.legend()
ax.set_xlabel(r'log $M_{200c}$ [km s$^{-1}$]')
ax.set_ylabel('redshift')
ax.set_ylim((0,zmax))
ax.set_zlabel(r'Data / Model')
ax.set_zlim((0.5,1.5))
#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir , "M200c-cumulative-function-allZ-modelRatio.pdf"))
p.clf()



sys.exit()
################################ Plot differential halo mass function  ################################


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_1Gpc"+".dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_2.5Gpc"+".dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),unpack=True)

#rhom_04 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_04])
#rhom_10 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_10])
#rhom_25 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_25])
#rhom_40 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_40])

rhom_04 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_04])
rhom_10 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_10])
rhom_25 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_25])
rhom_40 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_40])

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

#p.plot(n.log10(dat[0]),n.log10(dat[0]*dat[0]*dat[5]),'k--',lw=2)
s_04 = (z_04 <= 0.01)
y_04 = n.log10(yData_04*xData_04**2./rhom_04)
#y_04 = n.log10(yData_04*xData_04*xData_04)
p.plot(n.log10(xData_04[s_04]), y_04[s_04], marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)

s_10 = (z_10 == 0)
y_10 = n.log10(yData_10*xData_10**2./rhom_10)
#y_10 = n.log10(yData_10*xData_10*xData_10)
p.plot(n.log10(xData_10[s_10]),y_10[s_10], marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)

s_25 = (z_25 == 0)
y_25 = n.log10(yData_25*xData_25**2./rhom_25)
#y_25 = n.log10(yData_25*xData_25*xData_25)
p.plot(n.log10(xData_25[s_25]),y_25[s_25], marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)

s_40 = (z_40 == 0)
y_40 = n.log10(yData_40*xData_40**2./rhom_40)
#y_40 = n.log10(yData_40*xData_40*xData_40)
p.plot(n.log10(xData_40[s_40]),y_40[s_40], marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)

p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[(M^2/\rho_m) \; dn/dM]')
p.legend(loc=3)
p.grid()
p.savefig(join(Pdir , "M200c-diff-function-z0.pdf"))
p.clf()

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')

sc1 = ax.scatter(n.log10(xData_04),z_04,y_04, s=n.ones_like(z_04)*3, c='r', marker='o',label="SMD", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_10),z_10,y_10, s=n.ones_like(z_10)*3, c='c', marker='v',label="MDPL", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_25),z_25,y_25, s=n.ones_like(z_25)*3, c='m', marker='s',label="BigMD", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_40),z_40,y_40, s=n.ones_like(z_40)*3, c='b', marker='p',label="HMD", rasterized=True)
sc1.set_edgecolor('face')

ax.legend()
ax.set_xlabel(r'$log_{10}[M_{200c}/(h^{-1}M_\odot)]$')
ax.set_ylabel('redshift')
ax.set_ylim((0,zmax))
ax.set_zlabel(r'$log_{10}[(M^2/\rho_m) \; dn/dM]$')
ax.set_zlim((-4,0))
p.savefig(join(Pdir , "M200c-diff-function-allZ.pdf"))
p.clf()

