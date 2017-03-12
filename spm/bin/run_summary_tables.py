import os
from os.path import join
import numpy as n
import time

plates = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "catalogs", "plateNumberList"), unpack=True, dtype='str')

tts=[]
for plate in plates[10:]:
	t0 = time.time()
	os.system("python test.py "+plate+" stellarpop-m11-kroupa")
	tts.append(time.time()-t0)
	print tts[-1], "seconds"
	

plates = n.loadtxt( join(os.environ['EBOSSDR14_DIR'], "catalogs", "plateNumberList"), unpack=True, dtype='str')

tts=[]
for plate in plates[4:]:
	t0 = time.time()
	os.system("python test2.py "+plate+" stellarpop-m11-kroupa")
	tts.append(time.time()-t0)
	print tts[-1], "seconds"