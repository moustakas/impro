;+
; NAME:
;   CONVERT_TO_MAGGIES()
;
; PURPOSE:
;   Convert magnitudes stored in a data structure to "maggies" for use
;   with K-correct.
;
; INPUTS: 
;   data - data structure containing the observed photometry and
;          errors; for example, the structure might look like: [NGAL] 
;      .B - B-band magnitude
;      .V - V-band magnitude
;      .R - R-band magnitude
;      .B_err - corresponding uncertainty
;      .V_err - corresponding uncertainty
;      .R_err - corresponding uncertainty
;
;   bands - tag names in DATA corresponding to the observed photometry
;     (e.g., ['B','V','R']) [NBAND]
;
; OPTIONAL INPUTS: 
;   vega2ab - conversion from Vega to AB (default 0.0); note that
;     unless the photometry is already in AB, then you *must* give
;     VEGA2AB values [NBAND]
;   err_suffix - suffix to append to BANDS indicating the tag names of
;     the magnitude *errors* (default '_err')
;
; KEYWORD PARAMETERS: 
;   nanomaggies - convert the output 'maggies' and errors to
;     nanomaggies 
;
; OUTPUTS: 
;   maggies - maggies [NBAND,NGAL]
;
; OPTIONAL OUTPUTS:
;   maggies_invvar - inverse variance corresponding to MAGGIES [NBAND,NGAL]
;   mags           - input magnitudes [NBAND,NGAL]
;   err_mags       - corresponding uncertainties [NBAND,NGAL]
;   ivar_mags      - corresponding inverse variances [NBAND,NGAL]
;   badvalue       - value indicating a null measurement on input and
;     output; for some stupid reason this should be a negative value;
;     (default -999.0) 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Aug 04, U of A - written
;   jm06sep15nyu - added MAGS (AB) and ERR_MAGS outputs
;   jm08aug11nyu - documented
;
; Copyright (C) 2005, 2006, 2008, John Moustakas
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

function convert_to_maggies, data, bands, vega2ab, err_suffix=err_suffix, $
  maggies_invvar=maggies_invvar, mags=mags, err_mags=mags_err, $
  ivar_mags=mags_ivar, badvalue=badvalue, nanomaggies=nanomaggies
    
    ngalaxy = n_elements(data)
    if (ngalaxy eq 0L) then begin
       doc_library, 'convert_to_magggies'
       return, -1L
    endif

    nbands = n_elements(bands)
    if (n_elements(vega2ab) eq 0L) then begin
       vega2ab = fltarr(nbands)
    endif else begin
       if (nbands ne n_elements(vega2ab)) then begin
          splog, 'BANDS and VEGA2AB must match!'
          return, -1L
       endif
    endelse

    if (n_elements(err_suffix) eq 0L) then err_suffix = '_err'
    if (n_elements(badvalue) eq 0L) then badvalue = -999.0
    if (badvalue ge 0.0) then begin
       splog, 'BADVALUE must be less than zero!'
       return, -1L
    endif

    if keyword_set(nanomaggies) then $
      nanofactor = 1D9 else nanofactor = 1.0D
    
    mags = fltarr(nbands,ngalaxy)
    mags_err = mags
    mags_ivar = mags
    
    maggies = mags
    maggies_invvar = mags+1D-32
    
    for iband = 0L, nbands-1L do begin

       true = tag_exist(data,bands[iband],index=bandindx)
       true = tag_exist(data,bands[iband]+err_suffix,index=banderrindx)

;      bad = where((data.(bandindx) ge badvalue) or (data.(banderrindx) ge badvalue) or $
;        (data.(banderrindx) eq 0.0),nbad,comp=good,ncomp=ngood)
       good = where((data.(bandindx) gt badvalue) and (data.(banderrindx) gt badvalue) and $
         (data.(banderrindx) gt 0.0),ngood)
       if (ngood ne 0L) then begin
          mags[iband,good] = (data.(bandindx))[good] + vega2ab[iband] ; convert from Vega to AB
          mags_err[iband,good] = (data.(banderrindx))[good]
          mags_ivar[iband,good] = 1.0/mags_err[iband,good]^2.0
          maggies[iband,good] = nanofactor*10.0^(-0.4*mags[iband,good])
          maggies_invvar[iband,good] = 1.0/(0.4*alog(10.0)*$
            maggies[iband,good]*mags_err[iband,good])^2
       endif
       
    endfor

    zero = where(mags eq 0.0,nzero)
    if (nzero ne 0L) then begin
       mags[zero] = badvalue
       mags_err[zero] = badvalue
    endif
    
return, maggies
end

