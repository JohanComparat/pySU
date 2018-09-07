#!/bin/bash/
cd $PYSU_DIR/galaxy/bin_eBOSS_ELG

# stacks after the luminosity function from $LFMODELS_DIR
python create_stack_list_UV.py
# long operation
python stack_spectra_UV

# fast operation
python fit_UV_lines.py

python plot_UV_stack.py


