#!/bin/bash
#PBS -l walltime=200:00:00
#PBS -o casTab.o.$PBS_JOBID
#PBS -e casTab.e$PBS_JOBID
#PBS -M comparat@mpe.mpg.de

pscp -r comparat@login5.sciama.icg.port.ac.uk:/mnt/lustre/sdss-dr12/eBOSS-DR14/catalogs/*short.fits D:\data\SDSS_SPM\

pscp -r comparat@login5.sciama.icg.port.ac.uk:/mnt/lustre/sdss-dr12/SDSS-DR12/catalogs/*short.fits D:\data\SDSS_SPM\



cd $EBOSSDR14_DIR/catalogs

java -jar ~/stilts.jar tcat in=FireflyGalaxyKroupaEbossDR14.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaEbossDR14short.fits

java -jar ~/stilts.jar tcat in=FireflyGalaxySalpeterEbossDR14.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxySalpeterEbossDR14short.fits


java -jar ~/stilts.jar tcat in=FireflyGalaxyKroupaNoDustEbossDR14.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaNoDustEbossDR14short.fits

java -jar ~/stilts.jar tcat in=FireflyGalaxySalpeterNoDustEbossDR14.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxySalpeterNoDustEbossDR14short.fits
			
rm flyAllList
ls $EBOSSDR14_DIR/stellarpop-m11-kroupa-nodust/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-eboss14-kr-nd.fits
			
rm flyAllList
ls $EBOSSDR14_DIR/stellarpop-m11-salpeter/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-eboss14-ss.fits

rm flyAllList
ls $EBOSSDR14_DIR/stellarpop-m11-salpeter-nodust/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-eboss14-ss-nd.fits

	
cd $SDSSDR12_DIR/catalogs

java -jar ~/stilts.jar tcat in=FireflyGalaxyKroupaSdss26.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaSdss26short.fits
		
		
java -jar ~/stilts.jar tcat in=FireflyGalaxyKroupaNodustSdss26.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaNodustSdss26short.fits
		

java -jar ~/stilts.jar tcat in=FireflyGalaxySalpeterSdss26.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxySalpeterSdss26short.fits
		

java -jar ~/stilts.jar tcat in=FireflyGalaxySalpeterNodustSdss26.fits ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxySalpeterNodustSdss26short.fits
		
