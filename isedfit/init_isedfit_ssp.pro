function init_isedfit_ssp, npix=npix, nage=nage
; jm10jan22ucsd - initialize the SSP structure for the various
; BUILD_*_SSP routines    
    ssp = {$
      Z:     0.0,$
      age:   dblarr(nage),$
      mstar: fltarr(nage),$
      wave:  fltarr(npix),$
      flux:  fltarr(npix,nage)}
return, ssp
end
