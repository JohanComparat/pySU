"""
This script produces the stacks for emission line luminosity limited samples.
"""
import sys
from SpectraStacking import *

tobestacked = glob.glob( "/home/comparat/database/VVDS/products/emissionLineLuminosityFunctions/*/*.fits")

for el in tobestacked:
	print "--------------------------------------"
	print el, glob.glob(el[:-5]+"*stack*"),len(n.array(glob.glob(el[:-5]+"*stack*")))
	if el.find('stack') > 0 :
		continue
	elif len(n.array(glob.glob(el[:-5]+"*stack*"))) > 0:
		continue
	else :
		print el, "stacks !"
		st=SpectraStacking(el, Nspec = 100)
		st.stackSpectra()

tobestacked = glob.glob( "/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/*/*.fits")

for el in tobestacked:
	print "--------------------------------------"
	print el, glob.glob(el[:-5]+"*stack*"),len(n.array(glob.glob(el[:-5]+"*stack*")))
	if el.find('stack') > 0 :
		continue
	elif len(n.array(glob.glob(el[:-5]+"*stack*"))) > 0:
		continue
	else :
		print el, "stacks !"
		st=SpectraStacking(el, Nspec = 400)
		st.stackSpectra()

sys.exit()
tobestacked = glob.glob( "/home/comparat/database/VIPERS/products/emissionLineLuminosityFunctions/*/*.fits")

for el in tobestacked:
	if el.find('stack') > 0 :
		continue
	else :
		st=SpectraStacking(el)
		st.stackSpectra()

