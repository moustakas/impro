function peg2isedfit, peg, sfr_const=sfr_const
; convert the output from IM_READ_PEG() to an iSEDfit compatible
; format

    dist = 10.0*3.085678D18     ; fiducial distance [10 pc in cm]

    nage = n_elements(peg)
    nlines = peg[0].nlines
    npix = n_elements(peg[0].wave)

    ised = {$
      age:      fltarr(nage),$ ; [yr]
      wave:     fltarr(npix),$
      flux:     fltarr(npix,nage),$
      linewave: fltarr(nlines),$
      lineflux: fltarr(nlines,nage),$
      mgalaxy:  fltarr(nage),$
      mstar:    fltarr(nage),$
      mgas:     fltarr(nage),$
      zgas:     fltarr(nage),$
      sfr:      fltarr(nage),$
      nlyc:     fltarr(nage)}

    ised.age = peg.age*1E6
    ised.wave = peg[0].wave
; normalize by the *stellar* mass    
    for ii = 0L, nage-1 do begin
       ised.flux[*,ii] = peg[ii].flux/peg[ii].mstar/(4.0*!dpi*dist^2.0)         ; [erg/s/cm2/A/M_sun]
       ised.lineflux[*,ii] = peg[ii].lineflux/peg[ii].mstar/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/M_sun]
    endfor
    ised.linewave = peg[0].linewave
    ised.mgalaxy = peg.mgalaxy
    ised.mstar = peg.mstar
    ised.mgas = peg.mgas
    ised.zgas = peg.zgas
    ised.sfr = peg.sfr
    ised.nlyc = peg.nlyc

; scale by the continuous SFR back to 1 M_sun/yr, if necessary 
    if (n_elements(sfr_const) ne 0) then begin
       factor = 1.0/(sfr_const*1D-6)
       ised.flux = ised.flux*factor
       ised.lineflux = ised.lineflux*factor
       ised.sfr = ised.sfr*factor
; ...and probably more       
    endif
    
return, ised
end

