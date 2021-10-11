#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.

python convert_stack_2_qmost_template.py

cp /data43s/SDSS/stacks/X_AGN/*stitched*


"""

### import modules used below
import os
import numpy as np
import astropy
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
from os.path import join
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as n
# speclite
import speclite.filters
from speclite.filters import FilterConvolution
from speclite.filters import ab_reference_flux

### define the required columns for a spectral template
minimum_required_keywords = ["LAMBDA", "FLUX_DENSITY"]
allowed_units = [
    u.Unit("Angstrom", format="fits"),
    u.Unit("erg / (Angstrom cm2 s)", format="fits"),
]

output1_directory = join( os.environ['HOME'], "SDSS/stacks" )
input1  =  join(output1_directory,"X_AGN", "DR16_ELG-stitched-stack.fits")
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits")
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits")
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT1-DR16QSO-stitched-stack.fits")
input1 =  join(output1_directory,"X_AGN", "DR16LRG-stitched-stack.fits")

def write_qmost_template(input1):
    output1_png =  input1[:-5] + "-qmost-template.png"
    output1_pdf =  input1[:-5] + "-qmost-template.pdf"
    output1_tpl =  input1[:-5] + "-qmost-template.fits"

    print(f"\n*** the following files will be accessed with this tutorial:")
    print(f"\t- input: {input1}")
    print(f"\t- output: {output1_png}")
    print(f"\t- output: {output1_pdf}")

    hdu1 = fits.open(input1)
    print(f"\n*** loading FITS file: {input1}")
    print(fits.info(input1))
    if hdu1[0].header.get("EXTEND", False):
        # FITS tables should have the bintable as an extension
        found_axes1 = hdu1[1].columns.info(output=False)
        data1 = hdu1[1].data
    else:
        # this section
        found_axes1 = hdu1[0].columns.info(output=False)
        data1 = hdu1[0].data
    print("\n*** found the following axes:")
    for name,unit in zip(found_axes1['name'], found_axes1['unit']):
        print(f"\t{name} ({unit})")
        if not unit in allowed_units:
            print(f"\t\t***** WARNING: {unit} was not in the allowed units: {allowed_units}!")

    ### define bounds
    wavelength_min = 3000
    wavelength_max = 11000

    ### arrays below and above the data bounds
    wavelength_step = np.diff(data1['wavelength']).mean()
    print("min: ", data1['wavelength'].min())
    print("max: ", data1['wavelength'].max())
    print("step: ", wavelength_step)

    ### define a new table, based on the old one but with appended data
    if data1['wavelength'].max() < wavelength_max:
        start = data1['wavelength'].max() + wavelength_step
        stop = wavelength_max + wavelength_step
        wave_red = np.arange(start=start, stop=stop, step=wavelength_step)
        print("will append the following to the 'red' wavelength spectrum: ", wave_red)
        nrows1 = len(data1['wavelength'])
        nrows2 = len(wave_red)
        nrows = nrows1 + nrows2
        wavelength = n.hstack((data1['wavelength'], wave_red))
        flux_density = n.hstack((data1['medianStack'], n.ones_like(wave_red)*data1['medianStack'][-1]))
    else:
        wavelength = data1['wavelength']
        flux_density = data1['medianStack']

    # where the filters are
    path = os.path.join(
        os.environ['GIT_AGN_MOCK'],
        'data',
        'photometry',
        'filters')

    Tab_HSC_r = Table.read(os.path.join(path, 'hsc', 'r_HSC.txt'), format='ascii')
    ok = (
        Tab_HSC_r['col1'] > 5190) & (
            Tab_HSC_r['col1'] < 7300) & (
                Tab_HSC_r['col2'] > 0)
    Tab_HSC_r['col2'][(ok == False)] = 0.

    Pass_HSC_r = speclite.filters.FilterResponse(
        wavelength=Tab_HSC_r['col1'] * u.AA,
        response=Tab_HSC_r['col2'],
        meta=dict(
            group_name='HSC',
            band_name='r'))

    hsc_r_filter = speclite.filters.load_filters('HSC-r')

    # redshift=0.0
    ll = wavelength * u.AA
    flambda = flux_density * 1e-17 * u.erg * u.cm**(-2) * u.s**(-1) *  u.AA**(-1)
    # ll_rf = ll / (1 + redshift)
    ABmag_sdss = hsc_r_filter.get_ab_magnitudes(flambda, ll )['HSC-r'][0]
    ##
    def rescale_by(r_mag_out, redshift): return 10 ** ((r_mag_out + 48.6) / -2.5) / 10 ** ( (ABmag_sdss + 48.6) / -2.5)
    # rescaling values
    rsbs = rescale_by(14, 0.0)

    ### save to new file
    hdu_cols = fits.ColDefs([
        fits.Column(name="LAMBDA", format='D', unit='Angstrom', array = wavelength ),
        fits.Column(name="FLUX_DENSITY", format='D', unit='erg / (Angstrom * cm * cm * s)', array = rsbs * flux_density * 1e-17 )
        ])

    hdu = fits.BinTableHDU.from_columns(hdu_cols)

    hdu.name = 'SPECTRUM'
    hdu.header['MAG'] = 14

    outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])

    outf.writeto(output1_tpl, overwrite=True)
    print(output1_tpl, 'written')#

    fig1 = plt.figure(10)
    for subplot in [211, 212]:
        plt.subplot(subplot)
        label1 = os.path.basename(input1)[:-5]
        plt.plot(hdu.data['LAMBDA'], hdu.data['FLUX_DENSITY'], label=label1)
        if subplot == 211:
            plt.xlim(None, 9500)
            # add HRS overlay
            ys = plt.gca().get_ylim()
            plt.fill_betweenx(ys, x1=[3700, 3700], x2=[5000, 5000], alpha=0.1, color='b')
            plt.fill_betweenx(ys, x1=[5000, 5000], x2=[7000, 7000], alpha=0.1, color='g')
            plt.fill_betweenx(ys, x1=[7000, 7000], x2=[9500, 9500], alpha=0.1, color='r')
            plt.xscale('log')
            plt.yscale('log')
            # plt.fill_betweenx(ys, x1=[3926,3926],x2=[4355, 4355], alpha=0.1, color='b')
            # plt.fill_betweenx(ys, x1=[5160,5160],x2=[5730, 5730], alpha=0.1, color='g')
            # plt.fill_betweenx(ys, x1=[6100,6100],x2=[6790, 6790], alpha=0.1, color='r')
        else:
            plt.xlabel(f"{found_axes1['name'][0]} ({found_axes1['unit'][0]})")
            #plt.ylabel(f"{found_axes1['name'][1]} ({found_axes1['unit'][1]})")
            plt.ylabel(f"{found_axes1['name'][1]} (erg/s/cm**2/Angstrom)")  # manually modified, for shorter length
            plt.yscale('log')
            plt.xlim(3700, 9500)
    # create legend and show plot
    #plt.legend()
    plt.show()
    fig1.savefig(output1_png, bbox_inches='tight')
    # fig1.savefig(output1_pdf, bbox_inches='tight')
    plt.clf()

input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT2-highZ-stitched-stack.fits")
write_qmost_template(input1)
input1  =  join(output1_directory,"X_AGN", "DR16_ELG-stitched-stack.fits")
write_qmost_template(input1)
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits")
write_qmost_template(input1)
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits")
write_qmost_template(input1)
input1 =  join(output1_directory,"X_AGN", "ROSAT_AGNT1-DR16QSO-stitched-stack.fits")
write_qmost_template(input1)
input1 =  join(output1_directory,"X_AGN", "DR16LRG-stitched-stack.fits")
write_qmost_template(input1)

