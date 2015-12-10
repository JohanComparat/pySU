import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from mpl_toolkits.mplot3d import Axes3D

Pdir = "/Users/johancomparat/Documents/papers-reports/2015_LF_O2HbO3/figures/Products_Galaxies/emissionLineLuminosityFunctions"
# "/home/comparat/database/Products_Galaxies/emissionLineLuminosityFunctions/" # on eboss
lines = "H1_4862", "O3_5007", "O2_3728"

###################### O2 3728 #######################
###################### O2 3728 #######################
###################### O2 3728 #######################
line = lines[2]
files = n.array(glob.glob(join(Pdir,line, "*fitInfo.PKL")))
files.sort()
print files, line

lstar, phistar, alpha, sigma, redshift = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

lstar_err, phistar_err, alpha_err, sigma_err = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

for ii in range(len(files)):
    f=open(files[ii],'r')
    data = cPickle.load(f)
    f.close()
    zz = float(data[2].split('_')[-1][1:])/1000.
    if zz>0.5 :
        redshift[ii] = zz
    if zz<0.5 :
        redshift[ii] = 1.2
    if len(data)==3:
        lstar[ii], phistar[ii], alpha[ii], sigma[ii] = data[0]		
        lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, data[1][3][3]**0.5		
    if len(data)==4:
        lstar[ii], phistar[ii], alpha[ii] = data[0]
        sigma[ii] = data[3]		
        lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, 0.		
    if len(data)==5:
        lstar[ii], phistar[ii]  = data[0] 		
        sigma[ii], alpha[ii] = data[3], data[4]
        lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, 0., 0.		



Ns=3

pl1 = lambda x, a, b : a * x + b
plaw = lambda x, a, b : 10**(a * x + b)
pl2 = lambda x, a0, a1, a2 : a0+ a1 * x + a2*x*x

ok = (redshift>0.0)

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], lstar[ok], lstar_err[ok],fmt='none')
res, cov = curve_fit(pl1, redshift[ok], lstar[ok], sigma =  lstar_err[ok], p0 = (1,41) , maxfev = 5000000)
p.plot(redshift[ok], pl1(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\log_{10}(L_*)$')
p.xlabel('z')
p.grid()
p.ylim((40,42))
p.xlim((0.55,1.3))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "lstar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(L_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"


p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], phistar[ok], phistar_err[ok],fmt='none')
res, cov = curve_fit(plaw, redshift[ok], phistar[ok], sigma =  phistar_err[ok], p0 = (1,-3) , maxfev = 5000000)
p.plot(redshift[ok], plaw(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\Phi_*$')
p.xlabel('z')
p.grid()
p.yscale('log')
p.ylim((0.0008,0.02))
p.xlim((0.55,1.3))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "phistar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(\Phi_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"

print alpha, alpha_err
print sigma, sigma_err
sys.exit()

###################### O3 5007 #######################
###################### O3 5007 #######################
###################### O3 5007 #######################
line = lines[1]
files = n.array(glob.glob(join(Pdir,line, "*fitInfo.PKL")))
files.sort()
print files, line

lstar, phistar, alpha, sigma, redshift = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

lstar_err, phistar_err, alpha_err, sigma_err = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

for ii in range(len(files)):
	f=open(files[ii],'r')
	data = cPickle.load(f)
	f.close()
	redshift[ii] = float(data[2].split('_')[-1][1:])/1000.
	if len(data)==3:
		lstar[ii], phistar[ii], alpha[ii], sigma[ii] = data[0]		
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, data[1][3][3]**0.5		
	if len(data)==4:
		lstar[ii], phistar[ii], alpha[ii] = data[0]
		sigma[ii] = data[3]		
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, 0.		
	if len(data)==5:
		lstar[ii], phistar[ii]  = data[0] 		
		sigma[ii], alpha[ii] = data[3], data[4]
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, 0., 0.		



Ns=3

pl1 = lambda x, a, b : a * x + b
plaw = lambda x, a, b : 10**(a * x + b)
pl2 = lambda x, a0, a1, a2 : a0+ a1 * x + a2*x*x

ok = (redshift>0.0)

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], lstar[ok], lstar_err[ok],fmt='none')
res, cov = curve_fit(pl1, redshift[ok], lstar[ok], sigma =  lstar_err[ok], p0 = (1,41) , maxfev = 5000000)
p.plot(redshift[ok], pl1(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\log_{10}(L_*)$')
p.xlabel('z')
p.grid()
p.ylim((40,42))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "lstar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(L_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"


p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], phistar[ok], phistar_err[ok],fmt='none')
res, cov = curve_fit(plaw, redshift[ok], phistar[ok], sigma =  phistar_err[ok], p0 = (1,-3) , maxfev = 5000000)
p.plot(redshift[ok], plaw(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\Phi_*$')
p.xlabel('z')
p.grid()
p.yscale('log')
p.ylim((0.0008,0.02))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "phistar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(\Phi_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"

print alpha, alpha_err
print sigma, sigma_err
sys.exit()

################## H BETA ####################
################## H BETA ####################
################## H BETA ####################
################## H BETA ####################
line = lines[0]
files = n.array(glob.glob(join(Pdir,line, "*fitInfo.PKL")))
files.sort()

lstar, phistar, alpha, sigma, redshift = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

lstar_err, phistar_err, alpha_err, sigma_err = n.empty(len(files)), n.empty(len(files)), n.empty(len(files)), n.empty(len(files))

for ii in range(len(files)):
	f=open(files[ii],'r')
	data = cPickle.load(f)
	f.close()
	redshift[ii] = float(data[2].split('_')[-1][1:])/1000.
	if len(data)==3:
		lstar[ii], phistar[ii], alpha[ii], sigma[ii] = data[0]		
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, data[1][3][3]**0.5		
	if len(data)==4:
		lstar[ii], phistar[ii], alpha[ii] = data[0]
		sigma[ii] = data[3]		
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, data[1][2][2]**0.5, 0.		
	if len(data)==5:
		lstar[ii], phistar[ii]  = data[0] 		
		sigma[ii], alpha[ii] = data[3], data[4]
		lstar_err[ii], phistar_err[ii], alpha_err[ii], sigma_err[ii] = data[1][0][0]**0.5, data[1][1][1]**0.5, 0., 0.		



Ns=3

pl1 = lambda x, a, b : a * x + b
plaw = lambda x, a, b : 10**(a * x + b)
pl2 = lambda x, a0, a1, a2 : a0+ a1 * x + a2*x*x

ok = (redshift>0.3)

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], lstar[ok], lstar_err[ok])
res, cov = curve_fit(pl1, redshift[ok], lstar[ok], sigma =  lstar_err[ok], p0 = (1,41) , maxfev = 5000000)
p.plot(redshift[ok], pl1(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\log_{10}(L_*)$')
p.xlabel('z')
p.grid()
p.ylim((40,42))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "lstar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(L_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"


p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshift[ok], phistar[ok], phistar_err[ok])
res, cov = curve_fit(plaw, redshift[ok], phistar[ok], sigma =  phistar_err[ok], p0 = (1,-3) , maxfev = 5000000)
p.plot(redshift[ok], plaw(redshift[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel(r'$\Phi_*$')
p.xlabel('z')
p.grid()
p.yscale('log')
p.ylim((0.0008,0.02))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(join(Pdir ,line, "phistar-"+line+"-z-param-amplitude-evolution.pdf"))
p.clf()

print "log_{10}(\Phi_*) & =(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") \\\\"


sys.exit()

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[1][ok], diag_err[1][ok])
res, cov = curve_fit(pl1, redshifts[ok],results[1][ok], sigma =  diag_err[1][ok], p0 = (-0.5,3) , maxfev = 5000000)
p.plot(redshifts[ok], pl1(redshifts[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.ylabel(r'$V_0$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-v0-evolution.pdf")
p.clf()
print "V_0=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),")"



p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[2][ok], diag_err[2][ok])
res, cov = curve_fit(pl2, redshifts[ok],results[2][ok], sigma =  diag_err[2][ok], p0 = (-0.7,2,1) , maxfev = 5000000)
p.plot(redshifts[ok], pl2(redshifts[ok],res[0], res[1],res[2]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.legend()
p.ylabel(r'$\alpha$')
p.xlabel('z')
p.grid()
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2))+" "+ str(n.round(res[2],2)))#+" "+ str(n.round(res[3],2)) )
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-alpha-evolution.pdf")
p.clf()
print "\\alpha=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") z + (", n.round(res[2],Ns), r"\pm ", n.round(cov[2][2]**0.5,Ns),")z^2"

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[3][ok], diag_err[3][ok])
res, cov = curve_fit(pl1, redshifts[ok],results[3][ok], sigma =  diag_err[3][ok], p0 = (-0.5,3) , maxfev = 5000000)
p.plot(redshifts[ok], pl1(redshifts[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.legend()
p.ylabel(r'$\beta$')
p.xlabel('z')
p.grid()
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))#+" "+ str(n.round(res[2],2))+" "+ str(n.round(res[3],2)) )
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-beta-evolution.pdf")
p.clf()
print "\\beta=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") "


"""
id = n.argsort(redshifts)
for ii in id : #range(id): 2))#len(redshifts)):
    print n.round(redshifts[ii],Ns)," & ", n.round(results[0][ii],Ns), r"$\pm$ ", n.round(diag_err[0][ii],Ns)," & ", n.round(results[1][ii],Ns),r" $\pm$ ", n.round(diag_err[1][ii],Ns)," & ", n.round(results[2][ii],Ns),r"$\pm$ ", n.round(diag_err[2][ii],Ns)," & ", n.round(results[3][ii],Ns)," $\pm$ ", n.round(diag_err[3][ii],Ns)," \\\\"#, n.round(chi2Rs[ii],Ns)," \\\\"

"""