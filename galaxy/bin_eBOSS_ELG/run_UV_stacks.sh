#!/bin/bash/
cd $PYSU_DIR/galaxy/bin_eBOSS_ELG

# stacks after the luminosity function from $LFMODELS_DIR
python create_stack_list_UV.py
# long operation
python stack_spectra_UV.py

# fast operation
python fit_UV_lines.py

python plot_UV_stack.py



python create_stack_list_UV_AGN_t2.py
python stack_spectra_AGN.py

import os
import numpy as n
ll = os.listdir('.')
testN = n.array([ len(os.listdir(li)) for li in ll])
# 
# scp -r u0936736@eboss.sdss.org:/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_10_10/spectra/full/9376/*.fits 9376/
# 
# rsync -avze u0936736@eboss.sdss.org:/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_10_10/spectra/full/* .
# 
# (py36he1srv) (0)comparat@he1srv:~$ rsync -avze u0936736@eboss.sdss.org:/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_10_10/spectra/full/10* .
