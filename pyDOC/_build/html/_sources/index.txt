
Welcome to Skies and Universes python documentation
===================================================

This set of modules helps exploiting the database of galaxy spectra and Multidark N body simulation available at http://eboss.ft.uam.es/~comparat/ELG-db/database/, which contains spectroscopic data from the COSMOS, DEEP2, SDSS, VIPERS and VVDS surveys.

The folder containing the code is here http://eboss.ft.uam.es/~comparat/ELG-db/database/pySU/ . 
It contains two folders: pyGalaxy and pyMultidark. Please add these folders to your pythonpath to be able to use it.

Galaxy surveys
==============

Modules and classes to fit models to data
-----------------------------------------

.. toctree::
   :maxdepth: 2 

   LineFittingLibrary
   LineLuminosityFunction
   ModelSpectraStacks

Modules and classes to handle the data
--------------------------------------

.. toctree::
   :maxdepth: 2 

   GalaxySurveyDEEP2
   GalaxySpectrumDEEP2

   GalaxySurveyVIPERS
   GalaxySpectrumVIPERS

   GalaxySurveyVVDS
   GalaxySpectrumVVDS
   SpectraStacking

Support libraries
-----------------

.. toctree::
   :maxdepth: 2 

   lib_plot
   filterList
   lineListAir
   MiscellanousFunctionsLibrary

Scripts available
-----------------
 * calibrate_DEEP2_spectra.py : performs the flux calibration of DEEP2 spectra
 * fit_lines_DEEP2_fc_spectra.py : fits emission lines on the DEEP2 spectra
 * fit_lines_VIPERS_spectra.py : fits emission lines on the VIPERS spectra
 * fit_lines_VVDSWIDE_spectra.py : fits emission lines on the VVDS WIDE spectra
 * fit_lines_VVDSDEEP_spectra.py : fits emission lines on the VVDS DEEP spectra
 * fit_lines_VVDSUDEEP_spectra.py : fits emission lines on the VVDS UDEEP spectra
 * plotSurveys.py : does quality plots to check everything is fine in the survey

Luminosity functions
--------------------
Scripts that work around the luminosity function

 * estimate_line_LF.py : estimates the emission line luminosity function of the DEEP2, VVDS, and VIPERS surevy
 * stack_spectra.py : stacks the spectra issued from the luminosity function
 * fit_model_LF.py : fits models to the LFs
 * plotLFs.py : plots the LFs 

MultiDark N body simulations :
==============================

https://www.cosmosim.org/

Modules and classes :
---------------------

Mass function, lightcones, clustering ...

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

