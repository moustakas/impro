;+
; NAME:
;       TELLURIC_MASK()
;
; PURPOSE:
;       Mask telluric absorption features and flag wavelengths that
;       are uncontaminated by the absorption features.
;
; CALLING SEQUENCE:
;       mask = telluric_mask(wave,good=,bad=)
;
; INPUTS:
;       wave  - wavelength vector [Angstrom]
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       mask  - pixel mask for WAVE (1B is good and 0B is bad)  
;
; OPTIONAL OUTPUTS:
;       good      - indices in WAVE that have not been masked
;       bad       - indices in WAVE that have been masked
;
; COMMENTS:
;       See Figure 1 in Bessell (1999) for the strongest telluric
;       feature wavelengths.  Also see
;       http://www2.keck.hawaii.edu/inst/hires/makeewww/Atmosphere.
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 November 27, U of A
;       jm03dec24uofa, added [B,A,T]BAND and [B,A,T]FREEBAND keywords 
;       jm04sep10uofa - removed [B,A,T]BAND and [B,A,T]FREEBAND
;                       keywords; updated and generalized the telluric
;                       masking to match ICONSTRUCT_TELLURIC
;
; Copyright (C) 2003-2004, John Moustakas
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

function telluric_mask, wave, good=good, bad=bad

    nwave = n_elements(wave)
    if nwave eq 0L then begin
       print, 'Syntax - mask = telluric_mask(wave,good=,bad=)'
       return, -1
    endif

    mask = bytarr(nwave)

    tellbands1 = { TELLBAND, $
      twave1: 6850., $
      twave2: 6960., $
      cwave1: [6600., 6950., 0], $
      cwave2: [6860., 7200., 0] }
    tellbands2 = { TELLBAND, $
      twave1: 7150., $
      twave2: 7350., $
      cwave1: [7050., 7115., 7340.], $
      cwave2: [7160., 7130., 7440.] }
    tellbands3 = { TELLBAND, $
      twave1: 7560., $
      twave2: 7720., $
      cwave1: [7400., 7700., 0], $
      cwave2: [7580., 8000., 0] }
    tellbands4 = { TELLBAND, $
      twave1: 8105., $
      twave2: 8240., $
      cwave1: [8000., 8225., 0], $
      cwave2: [8105., 8325., 0] }
;   tellbands5 = { TELLBAND, $
;     twave1: 8530., $
;     twave2: 8865., $
;     cwave1: [8490., 8865., 0], $
;     cwave2: [8530., 8905., 0] }
;   tellbands6 = { TELLBAND, $
;     twave1: 8644., $
;     twave2: 8697., $
;     cwave1: [8604., 8697., 0], $
;     cwave2: [8644., 8737., 0] }

    tellbands = [tellbands1, tellbands2, tellbands3, tellbands4]

    for iband = 0L, n_elements(tellbands)-1L do begin

       for i = 0L, n_elements(tellbands[iband].twave1)-1L do $
         mask = mask or (wave ge tellbands[iband].twave1[i] and $
           wave le tellbands[iband].twave2[i])

    endfor

    mask = mask eq 0B ; 1B = good
    
    good = where(mask eq 1B,ngood,comp=bad,ncomp=nbad)
    
return, mask
end
