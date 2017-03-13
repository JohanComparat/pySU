import os
from os.path import join
import numpy as n
import time

plates = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "catalogs", "plateNumberList"), unpack=True, dtype='str')

for plate in plates:
	os.system("python create_summary_tables_sdss_dr12.py "+plate+" stellarpop-m11-kroupa-nodust")
	os.system("python create_summary_tables_sdss_dr12.py "+plate+" stellarpop-m11-salpeter")
	os.system("python create_summary_tables_sdss_dr12.py "+plate+" stellarpop-m11-salpeter-nodust")
	
sys.exit()

plates = n.loadtxt( join(os.environ['EBOSSDR14_DIR'], "catalogs", "plateNumberList"), unpack=True, dtype='str')

for plate in plates:
	os.system("python create_summary_tables_eboss_dr14.py "+plate+" stellarpop-m11-salpeter")
	os.system("python create_summary_tables_eboss_dr14.py "+plate+" stellarpop-m11-salpeter-nodust")
	os.system("python create_summary_tables_eboss_dr14.py "+plate+" stellarpop-m11-kroupa-nodust")
	
	
