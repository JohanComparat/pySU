#!/bin/bash
#PBS -l walltime=200:00:00
#PBS -o spallFly.o.$PBS_JOBID
#PBS -e spallFly.e$PBS_JOBID
#PBS -M comparat@mpe.mpg.de

cd $EBOSSDR14_DIR/catalogs
ls $EBOSSDR14_DIR/stellarpop-m11-kroupa/flyAll_catalogs/spFlyAll*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaEbossDR14.fits
			
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
ls $SDSSDR12_DIR/stellarpop-m11-kroupa/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
	ocmd='sort plate' \
	omode=out ofmt=fits out=FireflyGalaxyKroupaSdss26.fits
			
rm flyAllList
ls $SDSSDR12_DIR/stellarpop-m11-kroupa-nodust/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-sdss26-kr-nd.fits
			
rm flyAllList
ls $SDSSDR12_DIR/stellarpop-m11-salpeter/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-sdss26-ss.fits

rm flyAllList
ls $SDSSDR12_DIR/stellarpop-m11-salpeter-nodust/flyAll_catalogs/*.fits > flyAllList

java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
	ocmd='sort plate' \
	omode=out ofmt=fits out=spFlyAll-sdss26-ss-nd.fits

#cd $SDSSDR12_DIR/catalogs
#ls $SDSSDR12_DIR/stellarpop-m11-kroupa/flyAll_catalogs/*.fits > flyAllList

#java -jar ~/stilts.jar tcat in=@flyAllList ifmt=fits \
#	icmd='keepcols "specObjID mjd plate fiberID run1d run2d plug_ra plug_dec z z_err zwarning class subclass z_noqso z_Err_noqso zWarning_noqso class_noqso subClass_noqso age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus stellar_mass stellar_mass_err_plus	 stellar_mass_err_minus spm_EBV nComponentsSSP"' \
#	ocmd='sort plate' \
#	omode=out ofmt=fits out=spFlyAll-kr.fits
			
