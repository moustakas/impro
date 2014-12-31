;+
; NAME:
;   INIT_ISEDFIT_SSP()
; PURPOSE:
;   Initialize a structure for the various BUILD_*_SSP routines.
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jan 22, UCSD
;-

function init_isedfit_ssp, npix=npix, nage=nage
    ssp = {$
      Zmetal:         0.0,$
      age:   fltarr(nage),$
      mstar: fltarr(nage),$
      nlyc:  fltarr(nage),$
      wave:  fltarr(npix),$
      flux:  fltarr(npix,nage)}
return, ssp
end
