
# second step spectroscopic success rate

# opens parent plate and distributions for infinite SNR

# opens the plates one by one and produces the sample
platelistELG = [8763, 8765, 8767]
mjdlistELG = [57309, 57307, 57310]

ii=0
mjd = str(mjdlistELG[ii])
plate = str(plateListELG[ii])

hduSpPlate = fits.open(join(os.environ['EBOSS_ROOT'],"spectro","redux",os.environ['RUN2D'],plate,"spPlate-"+plate+"-"+mjd+".fits))
hduSpZELGflag = fits.open(join(os.environ['EBOSS_ROOT'],"spectro","redux",os.environ['RUN2D'],plate,os.environ['RUN1D'],"spZ_ELGflag-"+plate+"-"+mjd+".fits))
hduSpZline = fits.open(join(os.environ['EBOSS_ROOT'],"spectro","redux",os.environ['RUN2D'],plate,os.environ['RUN1D'],"spZline-"+plate+"-"+mjd+".fits))
hduSpZall = fits.open(join(os.environ['EBOSS_ROOT'],"spectro","redux",os.environ['RUN2D'],plate,os.environ['RUN1D'],"spZall-"+plate+"-"+mjd+".fits))
hduSpZbest = fits.open(join(os.environ['EBOSS_ROOT'],"spectro","redux",os.environ['RUN2D'],plate,os.environ['RUN1D'],"spZbest-"+plate+"-"+mjd+".fits))



sys.exit()

# simulations :

# complete set of targets (here ELGs) on the sky :
raMockAll, decMockAll,redshiftMockAll = n.loadtxt(join( os.environ['EBOSS_TARGET'], "elg", "mockCatalogs", "SHAM_norm-mean1e+12-sig500000000000.0-fsat0.225_ELGdensity_radecz.cat"), unpack=True)

# this is the same file as :
# input distribution of targets (here ELGs) on the sky :
# hduMock = fits.open(join( os.environ['EBOSS_TARGET'], "ebosstarget", "v0899", "ebossMock-v0001-elg.fits"))
# hduMock[1].data['RA']
# hduMock[1].data['DEC']
# hduMock[1].data['redshift']
# and also the same file as :
# hduMock = fits.open(join(os.environ['EBOSS_TARGET'], "ebosstarget", "v0899", "ebosstarget-v0899-elg.fits"))

# has been tiled under the chunk :
chunk = "eboss889"
tilingDir = join(os.environ['EBOSS_ROOT'],"ebosstilelist/trunk/outputs",chunk)
tilingMaskDir = join(tilingDir,"masks")

hduMockTiled = fits.open(join(tilingDir,"final-"+chunk+".fits"))
raMockTiled = hduMockTiled[1].data['RA']
decMockTiled = hduMockTiled[1].data['DEC']

maskGeometry = mangle.Mangle(join(tilingDir,"geometry-"+chunk+".ply"))

# compute TSR : N galaxy with fiber / N galaxy targeted in each polygon.
pID_mock = n.array([ maskGeometry.polyid( raMockAll[ii], decMockAll[ii]) for ii in range(len(raMockAll)) ])
pID_mock_tiled = n.array([ maskGeometry.polyid( raMockTiled[ii], decMockTiled[ii]) for ii in range(len(raMockTiled)) ])

TSR = n.empty(len(raMockAll))
for id in maskGeometry.polyids:
	Ntargeted = len((pID_mock == id).nonzero()[0])
	NwithFiber = len((pID_mock_tiled == id).nonzero()[0])
	TSR[(pID_mock == id)] = n.ones_like(TSR[(pID_mock == id)])* float(NwithFiber) /Ntargeted

TSR[(pID_mock == -1)] = n.zeros_like(TSR[(pID_mock == -1)]) 

raMockAll[(TSR==0)]

# construct randoms :

maskCenterPost = join(tilingMaskDir,"centerpost_maskboss889.fits")

