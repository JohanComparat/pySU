
"""
Library of function to create halo catalogs matched to a density.

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

"""
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d
import scipy.stats as st
import os
from os.path import join

# initialize data side

# what tracer we are dealing with :
tracer_dir = join(os.environ['DATA_DIR'],"ELG")

# where is the NZ :
nz_dir = join(tracer_dir, "observations/NZ")
NZfile = join(nz_dir, "nz-fisherGRIW1.dat")

# loads the NZ, needs to be per deg2
zminIN, zmaxIN, nGal_Deg2IN = n.loadtxt( NZfile, unpack = True, usecols = (0,1,2) )
ok = (nGal_Deg2IN>0) & (zmaxIN<1.25)
zmin, zmax, nGal_Deg2  =  zminIN[ok], zmaxIN[ok], nGal_Deg2IN[ok]

# where the outputs will be stored :
mockOutpuName = "mocks_fischerGRIW1"
mockOutput_dir = join(tracer_dir,mockOutpuName)
os.system('mkdir ' + mockOutput_dir)

print "outputs stored : ",mockOutput_dir
print "N(z) from ", NZfile


# initialize the lightcone to extract the mock from
lcDir = "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/lightcones/lc_square_0.1z1.4/"
lcFile = join(lcDir,"lightcone_MD_2.5Gpc_0.4z1.4.fits")

# loads the LC
hdu = fits.open(lcFile)
hdR = hdu[0].header
print lcFile, " loaded, columns:"
print hdu[1].data.dtype
# [('ra', '>f8'), ('dec', '>f8'), ('z_real_space', '>f8'), ('z_redshift_space', '>f8'), ('v_parallel', '>f8'), ('id', '>i8'), ('num_prog', '>i2'), ('pid', '>i8'), ('upid', '>i8'), ('mvir', '>f4'), ('rvir', '>f8'), ('rs', '>f8'), ('vrms', '>f8'), ('vmax', '>f8'), ('Spin', '>f4'), ('M200b', '>f8'), ('M200c', '>f8'), ('b_to_a', '>f8'), ('c_to_a', '>f8'), ('Halfmass_Radius', '>f4'), ('Macc', '>f4'), ('Mpeak', '>f4'), ('Vacc', '>f8'), ('Vpeak', '>f8'), ('Acc_Rate_Inst', '>f4'), ('Acc_Rate_100Myr', '>f4'), ('Acc_Rate_Mpeak', '>f4')]

# properties of the lightcone
area = (2*30.)**2.

# boolean arrays that discriminate central and sat halos
cen = (hdu[1].data['pid'] ==  -1)
sat = (cen ==  False)

#function to slide in redshift
slice_Z = lambda hdu,minz,maxz : (hdu[1].data['z_redshift_space']>= minz)&(hdu[1].data['z_redshift_space']<maxz)


def writerCats(name,idSel):
    """ writes the obtained mock catalog for clustering estimation. """
    print "writes ", name
    n.savetxt(mockOutput_dir+name+"_radecz.cat",n.transpose([hdu[1].data['ra'][idSel], hdu[1].data['dec'][idSel], hdu[1].data['z_redshift_space'][idSel]]),fmt = '%.8f %.8f %.5f')

def writerCatsAll(name,idSel):
    """ writes the obtained mock catalog and a catalog with all the columns"""
    print "writes ", name
    n.savetxt(mockOutput_dir+name+"_radecz.cat",n.transpose([hdu[1].data['ra'][idSel], hdu[1].data['dec'][idSel], hdu[1].data['z_redshift_space'][idSel]]))
    tbhdu = fits.BinTableHDU.from_columns(hdu[1].colmuns)
    tbhdu.data = tbhdu[idSel]
    prihdu = fits.PrimaryHDU(header=hdR)
    thdulist = fits.HDUList([prihdu, tbhdu])
    outPutFileName = mockOutput_dir+name+"_allCols.fits"
    os.system('rm '+outPutFileName)
    thdulist.writeto(outPutFileName)

def get_distrib_QTY(hdu, colN, zmin, zmax):
    IDh = n.arange(len(hdu[1].data[colN]))
    zsel = slice_Z(hdu,zmin,zmax)
    IDhz = IDh[zsel] # all ids in this redshift bin
    QTY = hdu[1].data[colN][zsel] # all QTY in this redshift bin
    nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
    p.clf()
    print zmin,zmax,len(IDhz)
    return IDhz,QTY,nn,bb

def get_distrib_QTY_cen(hdu, colN, zmin, zmax):
    IDh = n.arange(len(hdu[1].data[colN]))
    zsel = slice_Z(hdu,zmin,zmax)&(cen)
    IDhz = IDh[zsel] # all ids in this redshift bin
    QTY = hdu[1].data[colN][zsel] # all QTY in this redshift bin
    nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
    p.clf()
    print zmin,zmax,len(IDhz)
    return IDhz,QTY,nn,bb

def get_distrib_QTY_sat(hdu, colN, zmin, zmax):
    IDh = n.arange(len(hdu[1].data[colN]))
    zsel = slice_Z(hdu,zmin,zmax)&(sat)
    IDhz = IDh[zsel] # all ids in this redshift bin
    QTY = hdu[1].data[colN][zsel] # all QTY in this redshift bin
    nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
    p.clf()
    print zmin,zmax,len(IDhz)
    return IDhz,QTY,nn,bb

def sham(nGal,IDhz, QTY, nn,bb):
    mfc = interp1d(nn,(bb[1:]+bb[:-1])/2.)
    QTYmax = mfc(len(QTY))
    QTYmin = mfc(len(QTY)-nGal)
    qsel = (QTY>QTYmin)&(QTY<= QTYmax)
    IDhzq = IDhz[qsel]
    print zmin,zmax,nGal,len(IDhzq)
    return IDhzq

def shamIncomplete(incompFactor, nGal,IDhz, QTY, nn,bb):
    mfc = interp1d(nn,(bb[1:]+bb[:-1])/2.)
    mfcInv = interp1d((bb[1:]+bb[:-1])/2.,nn)
    QTYmaxAll = mfc(len(QTY))/incompFactor
    Nmax = mfcInv(QTYmaxAll)
    QTYmax = mfc(Nmax)
    QTYmin = mfc(Nmax-nGal)
    qsel = (QTY>QTYmin)&(QTY<= QTYmax)
    IDhzq = IDhz[qsel]
    return IDhzq

def sham_QTY_max(QTY_max, nGal,IDhz, QTY, nn,bb):
    mfc = interp1d(nn,(bb[1:]+bb[:-1])/2.)
    mfcInv = interp1d((bb[1:]+bb[:-1])/2.,nn)
    Nmax = mfcInv(QTY_max)
    QTYmax = mfc(Nmax)
    QTYmin = mfc(Nmax-nGal)
    qsel = (QTY>QTYmin)&(QTY<= QTYmax)
    IDhzq = IDhz[qsel]
    return IDhzq

def selectGaussian(position,scatter, nGal,IDhz, QTY, nn,bb):
    # constructs the QTY intervals around the distribution
    expected_cdf = lambda x : st.norm.cdf(x, loc = position, scale = scatter)
    interval  =  [ position - 9 * scatter , position + 9 * scatter]
    xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
    out = expected_cdf(xs)
    expected_cdf_inv = interp1d(out,xs)
    boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
    # gets the number of halos to select
    expected_cdf_tot = lambda x : nGal * st.norm.cdf(x, loc = position, scale = scatter)
    Up = expected_cdf_tot(boundaries[1:])
    Low = n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
    N2select = Up-Low
    print N2select,Up,Low
    # select in mass in the box
    qsels = n.array([ (QTY>boundaries[ii])&(QTY<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
    IDhzqAll = n.array([ IDhz[qs] for qs in qsels ])
    # random downsample to the N2select in each bin
    i = 0
    ids_selected = []
    for arr in IDhzqAll:
        random.shuffle(arr)
        ids_selected.append(arr[:N2select[i]])
        i+= 1

    ids_selected = n.hstack(( n.array(ids_selected) ))
    return ids_selected

def selectGaussian_fsat(position,scatter,fsat, nGal,IDhz_c, QTY_c,IDhz_s, QTY_s ):
    nSat = int(nGal*fsat)
    print nGal,nSat,fsat,position,scatter
    # constructs the QTY intervals around the distribution 
    expected_cdf = lambda x : st.norm.cdf(x, loc = position, scale = scatter)
    interval  =  [ position - 9 * scatter , position + 9 * scatter]
    xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
    out = expected_cdf(xs)
    expected_cdf_inv = interp1d(out,xs)
    boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
    # gets the number of halos to select the SAT
    print "sats"
    expected_cdf_s  =  lambda x : nSat * st.norm.cdf(x, loc = position, scale = scatter)
    Up_s  =  expected_cdf_s(boundaries[1:])
    Low_s  =  n.hstack(( 0., expected_cdf_s(boundaries[1:])[:-1] ))
    N2select_s  =  Up_s-Low_s
    # select in mass in the box
    qsels_s = n.array([ (QTY_s>boundaries[ii])&(QTY_s<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
    IDhzqAll_s = n.array([ IDhz_s[qs] for qs in qsels_s ])

    # random downsample to the N2select in each bin
    i = 0
    ids_selected_s = []
    for arr2 in IDhzqAll_s:
        random.shuffle(arr2)
        print len(arr2),int(N2select_s[i])
        ids_selected_s.append(arr2[:int(N2select_s[i])])
        i+= 1

    id_s = n.hstack((n.array(ids_selected_s)))

    nSatReal = len(id_s)
    nCen = nGal-nSatReal
    print nGal,nSat,nCen,fsat,position,scatter

    # gets the number of halos to select the CEN, compatible with the sat fraction to get the right density.
    print "centrals"
    expected_cdf_c  =  lambda x : nCen * st.norm.cdf(x, loc = position, scale = scatter)
    Up_c  =  expected_cdf_c(boundaries[1:])
    Low_c  =  n.hstack(( 0., expected_cdf_c(boundaries[1:])[:-1] ))
    N2select_c  =  Up_c-Low_c
    # select in mass in the box
    qsels_c = n.array([ (QTY_c>boundaries[ii])&(QTY_c<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
    IDhzqAll_c = n.array([ IDhz_c[qs] for qs in qsels_c ])

    # random downsample to the N2select in each bin
    i = 0
    ids_selected_c = []
    for arr in IDhzqAll_c:
        random.shuffle(arr)
        print len(arr),int(N2select_c[i])
        ids_selected_c.append(arr[:int(N2select_c[i])])
        i+= 1

    id_c = n.hstack((n.array(ids_selected_c)))

    print len(id_c),len(id_s)
    ids_selected = n.hstack((id_c,id_s ))
    print len(ids_selected)
    return ids_selected




def selectLogNorm(position,scatter, nGal,IDhz, QTY, nn,bb):
    # constructs the QTY intervals around the distribution
    expected_cdf = lambda x : st.lognorm.cdf(x, position, scatter)
    interval  =  [ position - 9 * scatter , position + 9 * scatter]
    xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
    out = expected_cdf(xs)
    expected_cdf_inv = interp1d(out,xs)
    boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
    # gets the number of halos to select
    expected_cdf_tot = lambda x : nGal * st.lognorm.cdf(x, position, scatter)
    Up = expected_cdf_tot(boundaries[1:])
    Low = n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
    N2select = Up-Low
    print N2select,Up,Low
    # select in mass in the box
    qsels = n.array([ (QTY>boundaries[ii])&(QTY<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
    IDhzqAll = n.array([ IDhz[qs] for qs in qsels ])
    # random downsample to the N2select in each bin
    i = 0
    ids_selected = []
    for arr in IDhzqAll:
        random.shuffle(arr)
        ids_selected.append(arr[:N2select[i]])
        i+= 1

    ids_selected = n.hstack(( n.array(ids_selected) ))
    return ids_selected

