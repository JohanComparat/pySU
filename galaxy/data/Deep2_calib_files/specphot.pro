function specphot, spec1d, hdr,rlambda,rresponse,ilambda,iresponse
;+
; NAME:
;       SPECPHOT
;
; PURPOSE:
;       To add f_nu and f_lambda to a spec1d structure.
;
;       SPECPHOT first interpolates over the 1d spectrum's bad
;       pixels, extends the red and blue ends of the spectrum, and
;       measures the total photon counts through the CFHT 12K R and I  
;       filters.  It uses those counts, combined with the photometry 
;       in the provided FITS header, to remap the spectrum from 
;       counts/hour to ergs/cm^2/s/Hz (for FLUXN) 
;       or ergs/cm^2/s/Angstrom (for FLUXL)     
;
; SYNTAX
;      calspec1d=specphot(spec1d, hdr,[rlambda, rresponse, ilambda, iresponse])
;
; CATEGORY:
;      spec1d
;
;
; CALLING SEQUENCE:
;
;      calspec1d=specphot(spec1d, hdr,[rlambda, rresponse, ilambda, iresponse])
;
; INPUTS:
; spec1d = a 1d spectrum structure, as read in by FILL_GAP. 
;          It is recommend that qe_correction2 be run on this structure
;          before SPECPHOT, in order to take out chip-to-chip response
;          variations and DEIMOS throughput. 
;  hdr = the FITS header corresponding to spec1d; can be read in via
;        FILL_GAP at the same time as spec1d using the HEADER keyword
;       
;
; OPTIONAL INPUTS:
; rlambda = an array containg the wavelengths for the CFHT 12K R
;           filter curve
; rresponse = an array containing the normalized response for the CFHT
;             12K K filter curve 
; ilambda = an array containg the wavelengths for the CFHT 12K I
;           filter curve
; iresponse = an array containing the the normalized response for the
;           CFHT 12K I filter curve 
;

; OUTPUTS:
; spec1d = a structure containing all the tags from the original
;          spec1d file, as well as: 
;          fluxn = a vector containing the flux of the spectrum in
;                    units of ergs cm^-2 s^-1 Hz^-1
;          fluxl = a vector containg the flux of the spectrum in 
;                   units of ergs cm^-2s^-1 Angstrom^-1
;          ivar_fluxn = a vector containing the inverse variance 
;                       for flux_n
;          ivar_fluxl = a vector containing the inverse variance 
;                       for flux_l
;          fluxflag   = 0 if the fluxes were well-defined; 
;                  1 if there were negative counts in R (in which 
;                     case I is used to calibrate); 
;                  2 for negative counts in I (so we used R only); 
;                  3 if negative in both (in which case no flux 
;                  was determined)      
;                  4 if this doesn't seem to be a full, 1200-line
;                  grating spectrum
; 
;   
;  PROCEDURES CALLED:
;        djs_maskinterp 
;
;   RESTRICTIONS
;
;      ONLY suitable for DEEP2 1200-line data with central wavelength ~7800 Angstroms 
;
; EXAMPLE:
;        calspec1d = specphot(spec1d,hdr,rlambda, rresponse, ilambda)
;
;
; MODIFICATION HISTORY:
;        APW 8/20/09
;-  


if (max(spec1d.lambda)-min(spec1d.lambda)) gt 3000. $
   OR mean(spec1d.lambda) lt 7300 OR mean(spec1d.lambda) gt 8300 $
    OR SXPAR(hdr,'MAGR') eq 0 OR SXPAR(hdr,'MAGI') eq 0 then begin
	print,'WARNING: SPECPHOT will only work for DEEP2-like data'
        print
	print,'Returning input array with FLUXFLAG = 4'
        RETURN, {SPEC:spec1d.spec, LAMBDA:spec1d.lambda, IVAR:spec1d.ivar,$
  CRMASK:spec1d.crmask,  BITMASK:spec1d.bitmask, ORMASK:spec1d.ormask,$
  INFOMASK:spec1d.infomask, NBADPIX:spec1d.nbadpix, OBJPOS:spec1d.objpos,$
  FWHM:spec1d.fwhm, NSIGMA:spec1d.nsigma, R1:spec1d.r1, R2:spec1d.r2, $
  SKYSPEC:spec1d.skyspec, IVARFUDGE:spec1d.ivarfudge,$ 
  FLUXL:0, IVAR_FLUXL:0, FLUXN:0, IVAR_FLUXN:0,FLUXFLAG:4}
endif


spec1d.spec = djs_maskinterp(spec1d.spec, dilate((spec1d.ivar eq 0) $
   OR (spec1d.spec eq 0),fltarr(5)+1),/const)


;denominates red and blue ends of the spectrum
nonzero = where(erode(spec1d.ivar gt 0,fltarr(21)+1))
lambdamin = min(spec1d.lambda[nonzero], minidx)

lambdamax = max(spec1d.lambda[nonzero], maxidx)



;gives the angstroms per pixel of the blue and red regions
pixelblue = median(deriv(spec1d.lambda[nonzero[minidx: minidx + 200]]))
nblue = (lambdamin - 5000) / pixelblue
nblue = round(nblue)

pixelred = median(deriv(spec1d.lambda[nonzero[maxidx - 200:maxidx]]))
nred = (10000 - lambdamax) / pixelred
nred = round(nred)

;gives the blue and red range of wavelengths 
bluelambda = (lambdamin - nblue*pixelblue) + findgen(nblue) * pixelblue 
redlambda = lambdamax + findgen(nred) * pixelred


;gives the median counts for the blue and red ranges
minblueidx = (where(spec1d.ivar gt 0 ))[ minidx]
medblue = median(spec1d.spec[nonzero[minidx+10:minidx + 200]])
maxredidx = (where(spec1d.ivar gt 0 ))[ maxidx]
medred = median(spec1d.spec[nonzero[maxidx - 200: maxidx]])



;stores an array of the wavelength ranges extending spectrum to the blueand red
templambda = [bluelambda, spec1d.lambda[minblueidx:maxredidx], redlambda]
;stores an array of the counts extending the spectrum to the blue and red
tempspec = [bluelambda*0 + medblue, $
    ivarsmooth(spec1d.spec[minblueidx:maxredidx], $
     spec1d.ivar[minblueidx:maxredidx],41),$
           redlambda*0 + medred] 


;checks the keyword settings, if empty, reads in photon counts and
;wavelengths through the CFHT 12K R and I filters
if n_elements(rlambda) eq 0 or n_elements(rresponse) eq 0 then $
  readcol, getenv('IDLSPEC1D_DIR')+'/etc/Rresponse.txt', f = 'D,D,D,D,D,D,D', Rlambda,a,b,c,d,e, Rresponse


if n_elements(ilambda) eq 0 or n_elements(iresponse) eq 0 then $
   readcol, getenv('IDLSPEC1D_DIR')+'/etc/Iresponse.txt' ,f = 'D,D,D,D,D,D,D', Ilambda,a,b,c,d,e, Iresponse 


;reinterpolate the probability of passing through the filter from the R and I filter curves
; onto the actual wavelength array of the spectrum
weightr=interpol(Rresponse,Rlambda, templambda)
weighti=interpol(iresponse,ilambda, templambda)


;calculates the normalized counts based on the probability
;distribution for each filter
countr = total(weightr * tempspec)/total(weightr)
counti = total(weighti * tempspec)/total(weighti)

;reads in the magnitude for each filter
magr = sxpar(hdr, 'MAGR')
magi = sxpar(hdr, 'MAGI')

;calculates the flux at that magnitude (in erg/s/cm^2/Hz)
fluxr = 10^((magr + 48.6d0)/( -2.5)) 
fluxi = 10^((magi + 48.6d0)/(-2.5))


;gives the flux per count
fpcr = fluxr / countr
fpci = fluxi / counti


;effective wavelength 
;effr = total(Ilambda*Rresponse) / total(Rresponse)
;effi = total(Rlambda*Iresponse) / total(Iresponse)
effr = 6599.0889
effi = 8135.4026


x = [effr, effi]
y = [fpcr, fpci]

;creates a logarithmic fit between the flux per count values at the
;effective wavelengths for each filter
params = linfit(alog(x),alog(y))
fluxn_corr = exp(params[0])* spec1d.lambda^params[1] 

if total(finite(params)) EQ 0 THEN BEGIN
   if fpcr lt 0 then fluxflag=1
   if fpci lt 0 then fluxflag=2
   if fpci lt 0 AND fpci lt 0 then fluxflag=3
   fluxn_corr= (fpcr > fpci) > 0
endif else fluxflag=0

;gives the flux distribution: flux/count * count
; fluxn is f_nu: ergs/s/cm^2/Hz
fluxn = fluxn_corr * spec1d.spec
ivar_fluxn=spec1d.ivar/fluxn_corr^2

;converts the fluxn to fluxl
; fluxl is f_lambda: ergs/s/cm^2/Angstrom
fluxl_corr = (3d18) / (spec1d.lambda)^2
fluxl=fluxn*fluxl_corr
ivar_fluxl=ivar_fluxn/fluxl_corr^2

if max(fluxn) eq 0 then ivar_fluxn=fluxn
if max(fluxl) eq 0 then ivar_fluxl=fluxl

; return spectrum with extra tags added
RETURN, {SPEC:spec1d.spec, LAMBDA:spec1d.lambda, IVAR:spec1d.ivar,$
  CRMASK:spec1d.crmask,  BITMASK:spec1d.bitmask, ORMASK:spec1d.ormask,$
  INFOMASK:spec1d.infomask, NBADPIX:spec1d.nbadpix, OBJPOS:spec1d.objpos,$
  FWHM:spec1d.fwhm, NSIGMA:spec1d.nsigma, R1:spec1d.r1, R2:spec1d.r2, $
  SKYSPEC:spec1d.skyspec, IVARFUDGE:spec1d.ivarfudge,$ 
  FLUXL:fluxl, IVAR_FLUXL:ivar_fluxl, FLUXN:fluxn, IVAR_FLUXN:ivar_fluxn,FLUXFLAG:fluxflag}


END
