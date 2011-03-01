;+
; NAME:
;	COMPUTE_REFRACTION()
;
; PURPOSE:
;	Compute the atmospheric refraction with respect to a reference
;	wavelength at a given airmass.
;
; CALLING SEQUENCE:
;
; INPUTS:
;	airmass - airmass of the observation
;
; OPTIONAL INPUTS:
;	refwave - reference wavelength [Angstrom; default 5000]
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	deltar  - predicted atmospheric refraction at AIRMASS relative
;                 to REFWAVE over the wavelength range [3000,10000]
;                 Angstrom
;
; OPTIONAL OUTPUTS:
;	lambda  - wavelength vector at which the refraction was
;                 computed
;
; COMMENTS:
;	Reference: Filippenko 1982, PASP, 94, 715
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Feb 18, U of A
;-

function compute_refraction, airmass, refwave=refwave, lambda=lambda

    nairmass = n_elements(airmass)
    if nairmass eq 0L then begin
       print, 'Syntax - deltar = compute_refraction(airmass,[refwave=,lambda=])'
       return, 0L
    endif

    if not keyword_set(refwave) then refwave = 5000.0

    if nairmass gt 1L then begin

       for i = 0L, nairmass-1L do begin

          dr = compute_refraction(airmass[i],refwave=refwave,lambda=lambda)
          if i eq 0L then deltar = dr else deltar = [ [deltar], [dr] ]

       endfor

       return, deltar
          
    endif
    
; output wavelength vector
    
    minlambda = 3000.0 ; [Angstrom]
    maxlambda = 1E4    ; 1 micron [Angstrom]
    dlambda = 50.0     ; Angstrom per pixel

    lambda = findgen((maxlambda-minlambda)/dlambda)*dlambda+minlambda
    
; compute the refraction as a function of wavelength
; -----------------------------------------------------------
    
; atmospheric conditions at an altitude of ~2 km:

    P = 600.0D ; pressure [mm Hg]
    T = 7.0D   ; temperature [degrees C]
    f = 8.0D   ; water vapor pressure [mm Hg]

; refractive index of dry air at sea level minus 1

    n_sea = 1D-6*((64.328D)+(29498.1D)/((146D)-(1D-4*lambda)^(-2.0))+$
                  (255.4D)/((41D)-(1D-4*lambda)^(-2.0)))

; correct for the lower ambient temperature and pressure at the
; altitude of a typical observatory (refractive index minus 1)

    n_TP = n_sea * (P*(1+((1.049D)-(0.0157D)*T)*1D-6*P) / $
                    ((720.883D)*((1D)+0.003661*T)))
    
; correct for water vapor in the atmosphere (refractive index minus 1)

    n = 1D-6 * (n_TP * 1D6 - ((f*(0.0624D)-(0.000680D)*((1D-4*lambda)^(-2.0))) / $
                              ((1D)+(0.003661D)*T)))
    
; compute the differential refraciton relative to REFWAVE at the given
; airmass
; -----------------------------------------------------------

    za = acos(1./airmass)*!radeg           ; zenith angle [degrees]
    n_refwave = interpol(n,lambda,refwave)
    
    deltar = 206265.0D * (n - n_refwave) * tan(za*!dtor)

return, deltar
end
