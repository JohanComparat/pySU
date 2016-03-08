import numpy as n
from scipy.interpolate import interp1d
import pyneb as pn

# Conversion from Morton (1991, ApJS, 77, 119) wavelength in Angstrom
# SDSS spectra are in the vacuum, therefore the ref wavelengths of the lines must be in the vacuum.
AIR = lambda VAC : VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2 + 2.76249e8 / VAC**4)
vacs=n.arange(1000,12000,0.01)
airs=AIR(vacs)
VAC = interp1d(airs,vacs)

# Wavelengths from pyNeb Atoms are in A in vacuum like the SDSS spectra. No need to convert.

C3 = pn.Atom('C', 3)
#C3.printIonic()
C3_1908=1/(C3.getEnergy(C3.getTransition(1908)[0])-C3.getEnergy(C3.getTransition(1908)[1]))

C4 = pn.Atom('C', 4)
#C4.printIonic()
C4_1548=1/(C4.getEnergy(C4.getTransition(1548)[0])-C4.getEnergy(C4.getTransition(1548)[1]))

O2 = pn.Atom('O', 2)
#O2.printIonic()
O2_3727=1/(O2.getEnergy(O2.getTransition(3727)[0])-O2.getEnergy(O2.getTransition(3727)[1]))
O2_3729=1/(O2.getEnergy(O2.getTransition(3729)[0])-O2.getEnergy(O2.getTransition(3729)[1]))
#O2=AIR((O2_3727+O2_3729)/2.)
O2_mean=(O2_3727*3.326568+O2_3729*3.324086)/(3.326568 + 3.324086)

Ne3 = pn.Atom('Ne',3)
#Ne3.printIonic()
Ne3_3869=1/(Ne3.getEnergy(Ne3.getTransition(3869)[0])-Ne3.getEnergy(Ne3.getTransition(3869)[1]))
Ne3_3968=1/(Ne3.getEnergy(Ne3.getTransition(3968)[0])-Ne3.getEnergy(Ne3.getTransition(3968)[1]))

O3 = pn.Atom('O', 3)
#O3.printIonic()
O3_4363=1/(O3.getEnergy(O3.getTransition(4363)[0])-O3.getEnergy(O3.getTransition(4363)[1]))
O3_4960=1/(O3.getEnergy(O3.getTransition(4960)[0])-O3.getEnergy(O3.getTransition(4960)[1]))
#O3_5007=AIR(1/(O3.getEnergy(O3.getTransition(5007)[0])-O3.getEnergy(O3.getTransition(5007)[1])))
O3_5007=1/(O3.getEnergy(O3.getTransition(5007)[0])-O3.getEnergy(O3.getTransition(5007)[1]))

O1 = pn.Atom('O', 1)
O1_5578=1/(O1.getEnergy(O1.getTransition(5578)[0])-O1.getEnergy(O1.getTransition(5578)[1]))
O1_6302=1/(O1.getEnergy(O1.getTransition(6302)[0])-O1.getEnergy(O1.getTransition(6302)[1]))
O1_6365=1/(O1.getEnergy(O1.getTransition(6365)[0])-O1.getEnergy(O1.getTransition(6365)[1]))

N2 = pn.Atom('N', 2)
#N2.printIonic()
N2_5756=1/(N2.getEnergy(N2.getTransition(5756)[0])-N2.getEnergy(N2.getTransition(5756)[1]))
N2_6549=1/(N2.getEnergy(N2.getTransition(6549)[0])-N2.getEnergy(N2.getTransition(6549)[1]))
N2_6585=1/(N2.getEnergy(N2.getTransition(6585)[0])-N2.getEnergy(N2.getTransition(6585)[1]))

S2 = pn.Atom('S', 2)
#S2.printIonic()
S2_6718=1/(S2.getEnergy(S2.getTransition(6718)[0])-S2.getEnergy(S2.getTransition(6718)[1]))
S2_6732=1/(S2.getEnergy(S2.getTransition(6732)[0])-S2.getEnergy(S2.getTransition(6732)[1]))

Ar3 = pn.Atom('Ar', 3)
#Ar3.printIonic()
Ar3_7137=1/(Ar3.getEnergy(Ar3.getTransition(7137)[0])-Ar3.getEnergy(Ar3.getTransition(7137)[1]))

# Wavelengths from pyNeb RecAtoms are in A in Air like the SDSS spectra. Conversion needed.

H1=pn.RecAtom('H',1) # Hydrogen Balmer series
H1_3970=VAC(H1.getWave(7,2))
H1_4102=VAC(H1.getWave(6,2))
H1_4341=VAC(H1.getWave(5,2))
H1_4862=VAC(H1.getWave(4,2))
#H1.getWave(4,2) #VAC(H1.getWave(4,2))
H1_6564=VAC(H1.getWave(3,2))

H1=pn.RecAtom('H',1) # Hydrogen Lyman series
H1_1216=VAC(H1.getWave(2,1))

He1=pn.RecAtom('He',1) # Helium

He2=pn.RecAtom('He',2) # Helium
He2_4686=VAC(He2.getWave(4,3))
He2_5411=VAC(He2.getWave(7,4))

# Limits for the 4000 A fit
dl4k=150
#intLim4k=n.array([3950-dl4k, 3950, 4050, 4050+dl4k])
intLim4k=n.array([3600-dl4k, 3600, 4140, 4140+dl4k])

# limits for th eUV luminosities fits
intLimUV=n.array([2000,2200,3000,3200,3400,3600,4100,4300,4500,4700])

# system at 2360
# cmin1,cmax1=2080.,2240.
em1=2326.7
abs1=2343.7
em2=2365.3

aTR=2370.

abs2=2374.3
abs3=2382.2
em3=2396.2
# cmin2,cmax2=2400.,2550.

a0s2360=n.array([em1,abs1,em2,abs2,abs3,em3])

# system at 2600
em1=2586.1
em2=2599.6

aTR=2606.

abs1=2612.5
abs2=2626.3
#cmin1,cmax1=2400.,2550.
#cmin2,cmax2=2650.,2770.
a0s2600=n.array([em1,em2,abs1,abs2])


# system at 2800
Mga=2795.
Mgb=2802.
aTR=2798.
#cmin1,cmax1=2650.,2770.
#cmin2,cmax2=2807., 2840.
a0s2800=n.array([Mga,Mgb])
 
# abs2852=3851.9
# cmin2,cmax2=2870.,3000.
