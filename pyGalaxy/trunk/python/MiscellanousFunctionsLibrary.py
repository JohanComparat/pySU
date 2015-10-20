"""
This script contains miscellaneous functions.
"""

from scipy.interpolate import interp1d
import numpy as n

kla=lambda ll :2.659 *(-2.156+1.509/ll-0.198/ll**2+0.011/ll**3 ) + 4.05
klb=lambda ll :2.659 *(-1.857+1.040/ll)+4.05

def CalzettiLaw(ll):
	"""
	Calzetti (2000) extinction law
	:param ll: wavelength in Angstrom
	"""
	if ll>6300:
		return klb(ll)
	if ll<=6300:
		return kla(ll)


dat = n.loadtxt("../data/calzettiLaw.txt",unpack=True)
calzettiLaw = interp1d(dat[0],dat[1])

