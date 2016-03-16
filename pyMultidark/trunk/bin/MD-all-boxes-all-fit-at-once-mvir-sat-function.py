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

mf = lambda v, A, v0, alpha, beta : 10**A * (v/10**v0)**beta * n.e**(- (v/10**v0)**alpha )

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
qty = "mvir"

f=open(join("outputs_mvir", "outputs-mvir-sat-fits-all-at-once.txt"),'w')

f.write( "we consider the "+type+" of "+qty+" of "+ cos+ "\n")
f.write( "in the redshift range"+str(zmin)+" "+str(zmax)+ "\n")
#f.write( "for the boxes"+dir_boxes+ "\n")
#f.write( zList_files
f.write( "within the following limits for each box"+str(qty_limits) +"\n")
f.write( "each box has a volume of"+str(volume_boxes)+ "Mpc3/h3" +"\n")

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Satellite-mvir_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Satellite-mvir_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Satellite-mvir_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Satellite-mvir_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= 0.0) & (z_04 <= zmax)
s_10 = (z_10 >= 0.0) & (z_10 <= zmax)
s_25 = (z_25 >= 0.0) & (z_25 <= zmax)
s_40 = (z_40 >= 0.0) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10] ))#, z_25[s_25], z_40[s_40]))
f.write( "all redshifts available:"+str(set(redshift))+"\n")
mvir = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10])) )#, xData_25[s_25], xData_40[s_40])))
f.write( "min and max masses available:"+str(n.min(mvir))+" "+str(n.max(mvir))+"\n")
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10])))#, yData_25[s_25], yData_40[s_40])))
f.write( "min and max Y available:"+str( n.min(yData))+" "+str( n.max(yData))+"\n")
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10])))#, yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
f.write( "min and max Y error available:"+str( n.min(yDataErr))+" "+str( n.max(yDataErr))+"\n")

#first guess

vcut0 = 12.753 # 13.8425
vcut1 = -1.0 # -0.8762
vcut2 =0.0 #  0.0414

b0 = -0.9 # -0.8749
b1 = 0.0 # -0.0124
b2 = 0.0 # -0.0090

a0 = 0.395# 0.5702
a1 = -0.03 # -0.072
a2 = 0.0 # 0.014

A0 = -3.78 # -4.0178
A1 = 1.0 # 0.8232
A2 = 0.0 # -0.0714


vfG = lambda v, z, A0, A1, vcut0, vcut1, a0, a1, b0, b1 : n.log10( 10**(A0 + A1 * z) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 +a1*z) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7])
chi2fun = lambda ps : n.sum( (vfGbis(mvir,redshift,ps) - yData)**2. / (yData*0.017)**2. )/(len(yData) - 8)
p1 = n.array([ A0, A1, vcut0, vcut1, a0, a1, b0, b1 ])

f.write( "looks for the optimum parameters"+"\n")
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

pOpt = res.x
cov = res.direc
chi2perpoint = lambda ps : (vfGbis(mvir,redshift,ps) - yData)**2. / (yData*0.017)**2. 
chi2pp = chi2perpoint(pOpt)

f.write( "A(z) & = "+str(n.round(pOpt[0],NDecimal))+ " + "+str(n.round(pOpt[1],NDecimal))+r'\times z \\'+"\n")
f.write( r" M_{cut}(z) & = "+str(n.round(pOpt[2],NDecimal))+"+ "+str(n.round(pOpt[3],NDecimal))+r'\times z \\'+"\n")
f.write( r" \alpha(z) & = "+str(n.round(pOpt[4],NDecimal))+" + "+str(n.round(pOpt[5],NDecimal))+r'\times z \\'+"\n") #+ '+str(a2)+r'\times z^2 \\'
f.write( r" \beta(z) & = "+str(n.round(pOpt[6],NDecimal))+" + "+str(n.round(pOpt[7],NDecimal))+r'\times z \\'+"\n")


xModel = n.arange(n.min(mvir),14,0.1)
X,Y = n.meshgrid(xModel,n.arange(0,n.max(redshift)+0.025,0.025))
Z = vfGbis(X,Y,pOpt)

n.savetxt(join("outputs_mvir","mvir-sat-cumulative-function-best_fit.txt"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )

#######################################################
# now plots the results of the fit
print "now plots the results of the fit"

vmax_mod, z_mod, n_mod = n.loadtxt(join("outputs_mvir","mvir-sat-cumulative-function-best_fit.txt"), unpack=True)

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
p.savefig(join("outputs_mvir","mvir-sat-cumulative-function-model.png"))
p.clf()


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(mvir, yData, c=redshift,s=5, marker='o',label="data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.ylim((-7,0))
p.xlim((9.5, 14))
p.grid()
p.savefig(join("outputs_mvir","mvir-sat-cumulative-function-data.png"))
p.clf()



p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(mvir, redshift, c=chi2pp,s=5, marker='o', rasterized=True, vmin = 0, vmax = 4)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("$\chi^2$")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r'z') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.xlim((9.5, 14))
p.grid()
p.savefig(join("outputs_mvir","mvir-sat-cumulative-function-chi2pp.png"))
p.clf()


