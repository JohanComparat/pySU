#!/bin/bash/

cd $PYSU_DIR/pyGalaxies/trunk/bin

./calibrate_DEEP2_spectra

./fit_lines_DEEP2_fc_spectra
./fit_lines_VIPERS_spectra

# from here :
./fit_lines_VVDSDEEP_spectra
./fit_lines_VVDSUDEEP_spectra
./fit_lines_VVDSWIDE_spectra

./compute_line_luminosities_DEEP2
./compute_line_luminosities_VIPERS
./compute_line_luminosities_VVDS
./plotSurveys

./estimate_line_LF
./fit_model_LF
./plotLFs

./stack_spectra
./fit_lines_stacks
./model_stacked_spectra
./model_stacked_spectra_VVDS
#./model_stacked_spectra_plotOnly

./plot_firefly_model

./compute_intrinsic_LFs
