;+
; NAME:
;   MAG2MAGGIES()
;
; PURPOSE:
;   Convert from maggies to AB magnitude.
;
; INPUTS: 
;   maggies - input maggies [NBAND,NOBJ]
;
; OPTIONAL INPUTS: 
;   ivarmaggies - corresponding inverse variance maggies [NBAND,NOBJ] 
;   nsigma - NSIGMA upper limit for measurements measured at <NSIGMA
;     significance (default 2)
;
; OUTPUTS: 
;   mag - output AB magnitudes [NBAND,NOBJ]
;
; OPTIONAL OUTPUTS:
;   magerr - approximate 1-sigma error (assumes the error is symmetric
;     in magnitude) [NBAND,NOBJ] 
;   lomagerr - "lower" magnitude error [NBAND,NOBJ] 
;   himagerr - "upper" magnitude error [NBAND,NOBJ] 
;   magnsigma - NSIGMA magnitude limit
;
; COMMENTS:
;   
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 May 24, NYU
;   jm12feb28ucsd - added upper limits and asymmetric error bars 
;
; Copyright (C) 2011, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function maggies2mag, maggies, ivarmaggies=ivarmaggies, magerr=magerr, $
  lomagerr=lomagerr, himagerr=himagerr, magnsigmag=magnsigma, nsigma=nsigma

    mag = maggies*0.0-99.0
    if n_elements(ivarmaggies) ne 0 then begin
       if n_elements(nsigma) eq 0 then nsigma = 2.0

       magerr = ivarmaggies*0.0-99.0
       lomagerr = ivarmaggies*0.0-99.0
       himagerr = ivarmaggies*0.0-99.0
       magnsigma = ivarmaggies*0.0-99.0

       snr = maggies*sqrt(ivarmaggies) ; fractional error
       good = where(snr ge nsigma,ngood)
       upper = where((ivarmaggies ne 0.0) and snr lt nsigma,nupper)
       nodata = where((ivarmaggies eq 0.0),nnodata)
;      upper = where((maggies ne 0.0) and (ivarmaggies ne 0.0) and snr lt nsigma,nupper)
;      nodata = where((maggies eq 0.0) and (ivarmaggies eq 0.0),nnodata)

; significant detections
       if (ngood ne 0L) then begin
          mag[good] = -2.5*alog10(maggies[good])

          fracerr = 1.0/snr[good]
          magerr[good] = 2.5*alog10(exp(1))*fracerr ; Symmetric in magnitude (approximation)
          lomagerr[good] = 2.5*alog10(1+fracerr)    ; Bright end (flux upper limit)
          himagerr[good] = 2.5*alog10(1-fracerr)    ; Faint end  (flux lower limit)
       endif
; NSIGMA upper limits
       if (nupper ne 0L) then begin
          magnsigma[upper] = +2.5*alog10(sqrt(ivarmaggies[upper])/nsigma) ; note "+" instead of 1/ivar
       endif
    endif else begin
       good = where(maggies gt 0.0,ngood)
       if (ngood ne 0L) then mag[good] = -2.5*alog10(maggies[good])
    endelse

;   if (ngood ne 0L) then magerr[good] = 1.0/(0.4*sqrt(ivarmaggies[good])*maggies[good]*alog(10.0))

return, mag
end
    
