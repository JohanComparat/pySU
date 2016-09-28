
Welcome to Skies and Universes python documentation
===================================================

This set of modules helps exploiting the database of galaxy spectra and Multidark N body simulation available at http://projects.ift.uam-csic.es/skies-universes/, which currently contains spectroscopic data from the DEEP2, VIPERS and VVDS surveys.
 
It contains two folders of interest: pyGalaxy and pyMultidark. 
Please add these folders to your pythonpath to be able to use it.

Galaxy surveys
==========

Modules and classes to handle the data
---------------------------------------------

.. toctree::
   :maxdepth: 2 

   GalaxySurveyDEEP2
   GalaxySpectrumDEEP2
   GalaxySurveyVIPERS
   GalaxySpectrumVIPERS
   GalaxySurveyVVDS
   GalaxySpectrumVVDS
   SpectraStacking

Modules and classes to fit models to data
-----------------------------------------------

.. toctree::
   :maxdepth: 2 

   LineFittingLibrary
   LineLuminosityFunction
   ModelSpectraStacks

Support libraries
-----------------

.. toctree::
   :maxdepth: 2 

   lib_plot
   filterList
   lineListAir
   MiscellanousFunctionsLibrary

Scripts to construct catalogs
---------------------------------

 * calibrate_DEEP2_spectra : performs the flux calibration of DEEP2 spectra
 * fit_lines_DEEP2_fc_spectra : fits emission lines on the DEEP2 spectra
 * fit_lines_VIPERS_spectra : fits emission lines on the VIPERS spectra
 * fit_lines_VVDSWIDE_spectra : fits emission lines on the VVDS WIDE spectra
 * fit_lines_VVDSDEEP_spectra : fits emission lines on the VVDS DEEP spectra
 * fit_lines_VVDSUDEEP_spectra : fits emission lines on the VVDS UDEEP spectra
 * compute_line_luminosities_DEEP2 : adds line luminosities in a given cosmology
 * compute_line_luminosities_VIPERS : adds line luminosities in a given cosmology + aperture correction
 * compute_line_luminosities_VVDS : idem
 * plotSurveys : produces summary plots related to the surveys

Scripts to measure luminosity functions
----------------------------------------------

 * estimate_line_LF : estimates the emission line luminosity function of the DEEP2, VVDS, and VIPERS surevy
 * fit_model_LF_single : fits models to a the data sets separately
 * fit_model_LF_all : fits themodels on all the data sets at once
 * plotLFs : plots fits and the LFs 

 Scripts to stack spectra
 ---------------------------
 
 * stack_spectra_DEEP2 or _VVDSDEEP : stacks the spectra issued from the luminosity function
 
MultiDark N body simulations :
===================

https://www.cosmosim.org/

Modules and classes :
---------------------

.. toctree::
   :maxdepth: 2 

   MultiDark
   LightconeCreation
   HaloSelection

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

