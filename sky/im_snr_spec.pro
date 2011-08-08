;+
; NAME:
;       IM_SNR_SPEC()
;
; PURPOSE:
;       Compute the exposure time to achieve a given S/N for
;       spectroscopic observations.
;
; CALLING SEQUENCE:
;       result = im_snr_spec(paramfile,f_obj,pars,_extra=extra,/silent)
;
; INPUTS:
;       paramfile - parameter file describing the characteristics of
;                   the detector, the wavelength of observation,
;                   etc. (see for example
;                   $IMPRO_DIR/etc/snr_params_template.txt) 
;       f_obj     - scalar or vector target flux density [erg/s/cm2/A] 
;
; OPTIONAL INPUTS:
;       lambda    - override the LAMBDA parameter (can be a vector)
;       extra     - override any parameters in PARAMFILE
;
; KEYWORD PARAMETERS:
;       silent    - suppress messages to STDOUT
;
; OUTPUTS:
;       result - output data structure
;          f_obj   - input F_OBJ (unmodified)
;          snr     - input desired S/N (unmodified)
;          time    - exposure time required to achieve desired S/N [minutes]
;          e_obj   - total number of object electrons in the aperture 
;          e_sky   - total number of sky electrons in the aperture 
;          e_dark  - total number of dark current electrons in the aperture 
;          e_noise - total number of read noise electrons in the aperture 
;
; OPTIONAL OUTPUTS:
;       pars - data structure of the parameters in PARAMFILE
;
; PROCEDURES USED:
;       READCOL, IM_READ_VEGA(), TABINV, SPLOG
;
; COMMENTS:
;       Remember, to convert from surface brightness to flux density,
;       subtract 2.5*alog10(omega) where omega is the area [arcsec2]
;       subtended by the aperture.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 05, U of A - written
;-

function im_snr_spec, paramfile, f_obj, pars, lambda=lambda, silent=silent, _extra=extra

    common snr_spec, vega
    
    nparamfile = n_elements(paramfile)
    nobj = n_elements(f_obj)
    nlambda = n_elements(lambda)

    if (nparamfile eq 0L) or (nobj eq 0L) then begin
       print, 'Syntax - im_snr_spec '
       return, -1L
    endif

    if (nlambda ne 0L) then if (nlambda ne nobj) then begin
       print, 'F_OBJ and LAMBDA must have the same number of elements.'
       return, -1L
    endif
    
; call this routine recursively    
    
    if (nobj gt 1L) then begin
       for iobj = 0L, nobj-1L do begin
          if (nlambda ne 0L) then lambda1 = lambda[iobj]
          result1 = im_snr_spec(paramfile,f_obj[iobj],pars,lambda=lambda1,silent=silent,_extra=extra)
          if (iobj eq 0L) then result = result1 else result = [ [result], [result1] ]
       endfor
       return, reform(result)
    endif

; parse the parameter file
    
    if (n_elements(pars) eq 0L) then begin

       if (not keyword_set(silent)) then splog, 'Parsing the parameter file.'
       if file_test(paramfile) eq 0L then begin
          splog, 'Parameter file '+paramfile+' not found.'
          return, -1L
       endif

       readcol, paramfile, name, value, format='A,D', /silent, $
         comment='#', delimiter=' '
       nname = n_elements(name)

       pars = {paramfile: paramfile}
       for i = 0L, nname-1L do pars = create_struct(pars,name[i],value[i])       

    endif

    if (nlambda ne 0L) then pars.lambda = lambda
    
    if (n_elements(extra) ne 0L) then begin
       
       if tag_exist(extra,'SNR') then pars.snr = extra.snr
       if tag_exist(extra,'LAMBDA') then pars.lambda = extra.lambda
       if tag_exist(extra,'LUNARAGE') then pars.lunarage = extra.lunarage
       if tag_exist(extra,'DIAMTER') then pars.diameter = extra.diameter
       if tag_exist(extra,'EPSILON') then pars.epsilon = extra.epsilon
       if tag_exist(extra,'ETA') then pars.eta = extra.eta
       if tag_exist(extra,'DLAMBDA') then pars.dlambda = extra.dlambda
       if tag_exist(extra,'RDNOISE') then pars.rdnoise = extra.rdnoise
       if tag_exist(extra,'DARKRATE') then pars.darkrate = extra.darkrate
       if tag_exist(extra,'PIXSIZE') then pars.pixsize = extra.pixsize
       if tag_exist(extra,'NBINY') then pars.nbiny = extra.nbiny
       if tag_exist(extra,'NBINX') then pars.nbinx = extra.nbinx
       if tag_exist(extra,'SLIT') then pars.slit = extra.slit
       if tag_exist(extra,'APERTURE') then pars.aperture = extra.aperture
       
    endif

; compute some convenient numbers

    area = !dpi*(pars.diameter/2.0)^2.0          ; collecting area [cm2]
    npix = pars.aperture/pars.pixsize/pars.nbiny ; number of object & sky spatial pixels [pixel]
    omega = pars.slit*pars.aperture              ; solid angle subtended by the object and sky [arcsec2]

    hplanck = 6.6260755D-27      ; planck constant [cgs]
    light = 2.9979250D18         ; speed of light [Angstrom/s]
    constant = 1.0/hplanck/light ; constant of the equations
    
; read the vega SED and interpolate at LAMBDA

    if (n_elements(vega) eq 0L) then vega = im_read_vega()
    f_vega = interpol(vega.flux,vega.wave,pars.lambda) ; [erg/s/cm2/A]

; compute the sky surface brightness given LAMBDA and LUNARAGE;
; adopted from http://www.astro.utoronto.ca/~patton/astro/mags.html  

    age = [0,3,7,10,14]               ; lunar age [days]
    weff = [3605,4415,5515,6605,7900] ; effective wavelengths [Angstrom]

    lunar = [$
      [22.0,21.5,19.9,18.5,17.0], $ ; U
      [22.7,22.4,21.6,20.7,19.5], $ ; B
      [21.8,21.7,21.4,20.7,20.0], $ ; V
      [20.9,20.8,20.6,20.3,19.9], $ ; R
      [19.9,19.9,19.7,19.5,19.2]  $ ; I
      ]

    tabinv, weff, pars.lambda, lindx
    tabinv, age, pars.lunarage, aindx
    sigma_sky = interpolate(lunar,aindx,lindx)

; sky signal [electon/s/pixel]
    
    m_sky = sigma_sky-2.5*alog10(omega) ; sky magnitude in the aperture [mag]
    f_sky = f_vega*10^(-0.4*m_sky)      ; sky flux density [erg/s/cm2/A]

    I_sky = constant*area*pars.epsilon*pars.eta*pars.lambda*(pars.nbinx*pars.dlambda)*f_sky

; object signal [electon/s/pixel]; solve the quadratic equation for
; the time, given the desired S/N

    if (f_obj lt 0.0) then begin

       splog, 'F_OBJ is negative!'
       I_obj = 0.0D
       time = 0.0D

    endif else begin
    
       I_obj = constant*area*pars.epsilon*pars.eta*pars.lambda*(pars.nbinx*pars.dlambda)*f_obj
    
       aquad = (I_obj/pars.snr)^2
       bquad = -(I_obj+I_sky+npix*pars.darkrate)
       cquad = -npix*pars.rdnoise^2
       
       time = (-bquad+sqrt(bquad^2-4*aquad*cquad))/(2*aquad)/60.0 ; [minutes]

    endelse

    result = {$
      f_obj:   f_obj, $
      snr:     pars.snr, $
      time:    time, $
      e_obj:   sqrt(I_obj*time), $
      e_sky:   sqrt(I_sky*time), $
      e_dark:  sqrt(npix*pars.darkrate*time), $
      e_noise: sqrt(npix)*pars.rdnoise}

    if (not keyword_set(silent)) then begin
       
       splog, 'To achieve a S/N of '+string(pars.snr,format='(I0)')+' requires '+$
         string(time,format='(G0.0)')+' minutes.'
       splog, 'Noise terms: '
       splog, '  Object      : '+string(result.e_obj,format='(G0.0)')+' electron'
       splog, '  Sky         : '+string(result.e_sky,format='(G0.0)')+' electron'
       splog, '  Dark Current: '+string(result.e_dark,format='(G0.0)')+' electron'
       splog, '  Read Noise  : '+string(result.e_noise,format='(G0.0)')+' electron'

    endif

return, result
end    
