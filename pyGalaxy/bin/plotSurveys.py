"""
This script produces quality plots to check that the Surveys and their catalogs are fine.
"""

from lib_plot import *
from lineListAir import *
SNlim = 5


from GalaxySurveyVVDS import *
plotDir="/home/comparat/database/VVDS/products/emissionLineLuminosityFunctions/plotSurvey/"

# UDEEP
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_UDEEP_summary.LFcatalog.fits")

zok=(survey.catalog['ZFLAGS']>1.9)&(survey.catalog['ZFLAGS']<9.1)

plotRaDecEBV(survey.catalog['ALPHA'][zok], survey.catalog['DELTA'][zok], survey.catalog['EBV_MW'][zok] , "vvdsudeep-ra-dec.pdf", plotDir)
plotZ_SSR(survey.catalog['Z'][zok],survey.catalog['SSR'][zok],"vvdsudeep-z-ssr.pdf",'spectroscopic success rate',plotDir)
plot_I_TSR(survey.catalog['MAGI'][zok],survey.catalog['TSR'][zok],"vvdsudeep-z-tsr.pdf",'target sampling rate',plotDir)

ok=(zok) & (survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O2_3728_flux'][ok], "vvdsudeep-f-O2_3728-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{3728}_{II}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_luminosity']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O2_3728_luminosity'][ok], "vvdsudeep-L-O2_3728-z.pdf",r'$log_{10}(L\, [O^{3728}_{II}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_EW']>0.)&(survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O2_3728_EW'][ok], "vvdsudeep-EW-O2_3728-z.pdf",r'$log_{10}(EW\,[O^{3728}_{II}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O3_5007_flux'][ok], "vvdsudeep-f-O3_5007-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{5007}_{III}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_luminosity']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O3_5007_luminosity'][ok], "vvdsudeep-L-O3_5007-z.pdf",r'$log_{10}(L\, [O^{5007}_{III}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_EW']>0.) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O3_5007_EW'][ok], "vvdsudeep-EW-O3_5007-z.pdf",r'$log_{10}(EW\,[O^{5007}_{III}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['H1_4862_flux'][ok], "vvdsudeep-f-H1_4862-z.pdf",r'$log_{10}(\mathrm{flux }\, H^{4861}_{\beta})$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_luminosity']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['H1_4862_luminosity'][ok], "vvdsudeep-L-H1_4862-z.pdf",r'$log_{10}(L\, H^{4861}_{\beta})$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_EW']>0.)& (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['H1_4862_EW'][ok], "vvdsudeep-EW-H1_4862-z.pdf",r'$log_{10}(EW\,H^{4861}_{\beta})$ [A]', plotDir)

# WIDE
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_WIDE_summary.LFcatalog.fits")

zok=(survey.catalog['ZFLAGS']>1.9)&(survey.catalog['ZFLAGS']<9.1)
sel1=(survey.catalog['alpha']<170)&(zok)
sel2=(survey.catalog['delta']>3)&(zok)
sel3=(survey.catalog['alpha']>300)&(zok)
sel=n.array([sel1,sel2,sel3])

plotZ_SSR(survey.catalog['Z'][zok],survey.catalog['SSR'][zok],"vvdswide-z-ssr.pdf",'spectroscopic success rate',plotDir)
plot_I_TSR(survey.catalog['MAGI'][zok],survey.catalog['TSR'][zok],"vvdswide-z-tsr.pdf",'target sampling rate',plotDir)

for i in range(3):
	plotRaDecEBV(survey.catalog['ALPHA'][sel[i]], survey.catalog['DELTA'][sel[i]], survey.catalog['EBV_MW'][sel[i]] , "vvdswide-f"+str(i)+"-ra-dec.pdf", plotDir)

ok=(zok) & (survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O2_3728_flux'][ok], "vvdswide-f-O2_3728-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{3728}_{II}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_luminosity']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O2_3728_luminosity'][ok], "vvdswide-L-O2_3728-z.pdf",r'$log_{10}(L\, [O^{3728}_{II}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_EW']>0.)&(survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O2_3728_EW'][ok], "vvdswide-EW-O2_3728-z.pdf",r'$log_{10}(EW\,[O^{3728}_{II}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O3_5007_flux'][ok], "vvdswide-f-O3_5007-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{5007}_{III}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_luminosity']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O3_5007_luminosity'][ok], "vvdswide-L-O3_5007-z.pdf",r'$log_{10}(L\, [O^{5007}_{III}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_EW']>0.) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O3_5007_EW'][ok], "vvdswide-EW-O3_5007-z.pdf",r'$log_{10}(EW\,[O^{5007}_{III}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['H1_4862_flux'][ok], "vvdswide-f-H1_4862-z.pdf",r'$log_{10}(\mathrm{flux }\, H^{4861}_{\beta})$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_luminosity']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['H1_4862_luminosity'][ok], "vvdswide-L-H1_4862-z.pdf",r'$log_{10}(L\, H^{4861}_{\beta})$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_EW']>0.)& (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['H1_4862_EW'][ok], "vvdswide-EW-H1_4862-z.pdf",r'$log_{10}(EW\,H^{4861}_{\beta})$ [A]', plotDir)

# DEEP
survey = GalaxySurveyVVDS(redshift_catalog="VVDS_DEEP_summary.LFcatalog.fits")

zok=(survey.catalog['ZFLAGS']>1.9)&(survey.catalog['ZFLAGS']<9.1)

plotRaDecEBV(survey.catalog['ALPHA'][zok], survey.catalog['DELTA'][zok], survey.catalog['EBV_MW'][zok] , "vvdsdeep-ra-dec.pdf", plotDir)
plotZ_SSR(survey.catalog['Z'][zok],survey.catalog['SSR'][zok],"vvdsdeep-z-ssr.pdf",'spectroscopic success rate',plotDir)
plot_I_TSR(survey.catalog['MAGI'][zok],survey.catalog['TSR'][zok],"vvdsdeep-z-tsr.pdf",'target sampling rate',plotDir)

ok=(zok) & (survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O2_3728_flux'][ok], "vvdsdeep-f-O2_3728-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{3728}_{II}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_luminosity']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O2_3728_luminosity'][ok], "vvdsdeep-L-O2_3728-z.pdf",r'$log_{10}(L\, [O^{3728}_{II}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_EW']>0.)&(survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O2_3728_EW'][ok], "vvdsdeep-EW-O2_3728-z.pdf",r'$log_{10}(EW\,[O^{3728}_{II}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['O3_5007_flux'][ok], "vvdsdeep-f-O3_5007-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{5007}_{III}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_luminosity']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['O3_5007_luminosity'][ok], "vvdsdeep-L-O3_5007-z.pdf",r'$log_{10}(L\, [O^{5007}_{III}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_EW']>0.) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['O3_5007_EW'][ok], "vvdsdeep-EW-O3_5007-z.pdf",r'$log_{10}(EW\,[O^{5007}_{III}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Flux(survey.catalog['Z'][ok],survey.catalog['H1_4862_flux'][ok], "vvdsdeep-f-H1_4862-z.pdf",r'$log_{10}(\mathrm{flux }\, H^{4861}_{\beta})$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_luminosity']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['Z'][ok],survey.catalog['H1_4862_luminosity'][ok], "vvdsdeep-L-H1_4862-z.pdf",r'$log_{10}(L\, H^{4861}_{\beta})$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_EW']>0.)& (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_EW(survey.catalog['Z'][ok],survey.catalog['H1_4862_EW'][ok], "vvdsdeep-EW-H1_4862-z.pdf",r'$log_{10}(EW\,H^{4861}_{\beta})$ [A]', plotDir)


from GalaxySurveyDEEP2 import *
plotDir="/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/plotSurvey/"
survey = GalaxySurveyDEEP2(redshift_catalog="zcat.deep2.dr4.v2.LFcatalog.fits", calibration = False)

zok=(survey.catalog['ZQUALITY']>0.9)
sel1=(survey.catalog['DEC']>50)&(zok)
sel2=(survey.catalog['DEC']<50)&(survey.catalog['DEC']>30)&(zok)
sel3=(survey.catalog['DEC']<30)&(survey.catalog['RA']>200)&(zok)
sel4=(survey.catalog['DEC']<30)&(survey.catalog['RA']<200)&(zok)
sel=n.array([sel1,sel2,sel3,sel4])

plotZ_SSR(survey.catalog['ZBEST'][zok],survey.catalog['SSR'][zok],"deep2-z-ssr.pdf",'spectroscopic success rate',plotDir)
plot_I_TSR(survey.catalog['MAGR'][zok],survey.catalog['TSR'][zok],"deep2-z-tsr.pdf",'target sampling rate',plotDir)

for i in range(4):
	plotRaDecEBV(survey.catalog['RA'][sel[i]], survey.catalog['DEC'][sel[i]], survey.catalog['SFD_EBV'][sel[i]] , "deep2-f"+str(i)+"-ra-dec.pdf", plotDir)


ok=(zok) & (survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Flux(survey.catalog['ZBEST'][ok],survey.catalog['O2_3728_flux'][ok], "deep2-f-O2_3728-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{3728}_{II}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_luminosity']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['ZBEST'][ok],survey.catalog['O2_3728_luminosity'][ok], "deep2-L-O2_3728-z.pdf",r'$log_{10}(L\, [O^{3728}_{II}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_EW']>0.)&(survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_EW(survey.catalog['ZBEST'][ok],survey.catalog['O2_3728_EW'][ok], "deep2-EW-O2_3728-z.pdf",r'$log_{10}(EW\,[O^{3728}_{II}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Flux(survey.catalog['ZBEST'][ok],survey.catalog['O3_5007_flux'][ok], "deep2-f-O3_5007-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{5007}_{III}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_luminosity']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['ZBEST'][ok],survey.catalog['O3_5007_luminosity'][ok], "deep2-L-O3_5007-z.pdf",r'$log_{10}(L\, [O^{5007}_{III}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_EW']>0.) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_EW(survey.catalog['ZBEST'][ok],survey.catalog['O3_5007_EW'][ok], "deep2-EW-O3_5007-z.pdf",r'$log_{10}(EW\,[O^{5007}_{III}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Flux(survey.catalog['ZBEST'][ok],survey.catalog['H1_4862_flux'][ok], "deep2-f-H1_4862-z.pdf",r'$log_{10}(\mathrm{flux }\, H^{4861}_{\beta})$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_luminosity']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['ZBEST'][ok],survey.catalog['H1_4862_luminosity'][ok], "deep2-L-H1_4862-z.pdf",r'$log_{10}(L\, H^{4861}_{\beta})$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_EW']>0.)& (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_EW(survey.catalog['ZBEST'][ok],survey.catalog['H1_4862_EW'][ok], "deep2-EW-H1_4862-z.pdf",r'$log_{10}(EW\,H^{4861}_{\beta})$ [A]', plotDir)


from GalaxySurveyVIPERS import *
plotDir="/home/comparat/database/VIPERS/products/emissionLineLuminosityFunctions/plotSurvey/"

# VIPERS
survey = GalaxySurveyVIPERS(redshift_catalog="VIPERS_W14_summary_v1.LFcatalog.fits")

zok=(survey.catalog['zflg']>0.9)&(survey.catalog['zflg']<100)

plotRaDecEBV(survey.catalog['ALPHA'][zok], survey.catalog['DELTA'][zok], survey.catalog['EBV_MW'][zok] , "vipers-ra-dec.pdf", plotDir)
plotZ_SSR(survey.catalog['zspec'][zok],survey.catalog['SSR'][zok],"vipers-z-ssr.pdf",'spectroscopic success rate',plotDir)
plot_I_TSR(survey.catalog['MAGI'][zok],survey.catalog['TSR'][zok],"vipers-z-tsr.pdf",'target sampling rate',plotDir)

ok=(zok) & (survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Flux(survey.catalog['zspec'][ok],survey.catalog['O2_3728_flux'][ok], "vipers-f-O2_3728-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{3728}_{II}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_luminosity']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['zspec'][ok],survey.catalog['O2_3728_luminosity'][ok], "vipers-L-O2_3728-z.pdf",r'$log_{10}(L\, [O^{3728}_{II}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O2_3728_EW']>0.)&(survey.catalog['O2_3728_flux']>0.)& (survey.catalog['O2_3728_flux'] > SNlim*survey.catalog['O2_3728_fluxErr']) & (survey.catalog['O2_3728_fluxErr']>0.)
plotZ_EW(survey.catalog['zspec'][ok],survey.catalog['O2_3728_EW'][ok], "vipers-EW-O2_3728-z.pdf",r'$log_{10}(EW\,[O^{3728}_{II}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Flux(survey.catalog['zspec'][ok],survey.catalog['O3_5007_flux'][ok], "vipers-f-O3_5007-z.pdf",r'$log_{10}(\mathrm{flux }\, [O^{5007}_{III}])$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_luminosity']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['zspec'][ok],survey.catalog['O3_5007_luminosity'][ok], "vipers-L-O3_5007-z.pdf",r'$log_{10}(L\, [O^{5007}_{III}])$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['O3_5007_EW']>0.) & (survey.catalog['O3_5007_flux']>0.)& (survey.catalog['O3_5007_flux'] > SNlim*survey.catalog['O3_5007_fluxErr']) & (survey.catalog['O3_5007_fluxErr']>0.)
plotZ_EW(survey.catalog['zspec'][ok],survey.catalog['O3_5007_EW'][ok], "vipers-EW-O3_5007-z.pdf",r'$log_{10}(EW\,[O^{5007}_{III}])$ [A]', plotDir)

ok=(zok) & (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Flux(survey.catalog['zspec'][ok],survey.catalog['H1_4862_flux'][ok], "vipers-f-H1_4862-z.pdf",r'$log_{10}(\mathrm{flux }\, H^{4861}_{\beta})$ [erg s$^{-1}$ cm$^{-2}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_luminosity']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_Luminosity(survey.catalog['zspec'][ok],survey.catalog['H1_4862_luminosity'][ok], "vipers-L-H1_4862-z.pdf",r'$log_{10}(L\, H^{4861}_{\beta})$ [erg s$^{-1}$]', plotDir)
ok=(zok) & (survey.catalog['H1_4862_EW']>0.)& (survey.catalog['H1_4862_flux']>0.)& (survey.catalog['H1_4862_flux'] > SNlim*survey.catalog['H1_4862_fluxErr']) & (survey.catalog['H1_4862_fluxErr']>0.)
plotZ_EW(survey.catalog['zspec'][ok],survey.catalog['H1_4862_EW'][ok], "vipers-EW-H1_4862-z.pdf",r'$log_{10}(EW\,H^{4861}_{\beta})$ [A]', plotDir)
