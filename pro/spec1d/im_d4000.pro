;+
; NAME:
;   IM_D4000()
;
; PURPOSE:
;   Compute the Balogh+99 definition of the 4000-A break. 
;
; INPUTS: 
;   wave - input *rest-frame* wavelength vector [NPIX, Angstroms]
;   flux - corresponding flux density [NPIX, erg/s/cm/A]
;
; OPTIONAL INPUTS: 
;   ivar - inverse variance spectrum [NPIX, 1/(erg/s/cm/A)^2] (set
;     IVAR=0 for crummy pixels) 
;
; KEYWORD PARAMETERS: 
;   bruzual83 - compute the Bruzual (1983) definition of the 4000-A
;     break
;   debug - make a debugging plot and wait for a keystroke
;
; OUTPUTS: 
;   d4000 - the 4000-A break strength (dimensionless)
;
; OPTIONAL OUTPUTS:
;   d4000_err - uncertainty on D4000 (if IVAR is given, otherwise
;     -2.0) 
;
; COMMENTS:
;   Basically wraps on IABSLINEEW()
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 27, NYU
;   jm15apr20siena - changed FERR optional input to IVAR to conform to
;     the corresponding change in IABSLINEEW()
;
; Copyright (C) 2009, 2012, John Moustakas
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

function im_d4000, wave, flux, ivar=ivar, d4000_err=d4000_err, $
  bruzual83=bruzual83, debug=debug

    npix = n_elements(wave)
    if (npix eq 0) then begin
       doc_library, 'im_d4000'
       return, -1
    endif

    if (npix ne n_elements(flux)) then begin
       splog, 'Dimensions of WAVE and FLUX do not match!'
       return, -1
    endif
    
    if keyword_set(bruzual83) then begin
       wbvec = [3750.0,3950.0]
       wrvec = [4050.0,4250.0]
    endif else begin
       wbvec = [3850.0,3950.0]
       wrvec = [4000.0,4100.0]
    endelse

    llimit = total(wbvec)/2.0
    ulimit = total(wrvec)/2.0
    lwidth = wbvec[1]-wbvec[0]
    uwidth = wrvec[1]-wrvec[0]
    midwave = total([llimit,ulimit])/2.0 ; central wavelength

    cbreak = iabslineew(wave,flux,midwave,ivar=ivar,llimit=llimit,$
      lwidth=lwidth,ulimit=ulimit,uwidth=uwidth,label=blabel,/noline,$
      absplot=cbreakplot,debug=debug,silent=silent,_extra=extra)

    if (cbreak.cratio gt 0.0) then begin
       d4000 = cbreak.cratio
       d4000_err = cbreak.cratio_err
    endif else begin
       d4000 = -1.0
       d4000_err = -1.0
    endelse

    if (n_elements(ivar) eq 0) then d4000_err = -2.0    
    
return, d4000
end
    
