;+
; NAME:
;   CFHTLS_TO_MAGGIES
;
; PURPOSE:
;   Convert a Gwyn-style CFHTLS photometric catalog to AB maggies. 
;
; INPUTS: 
;   cat - input photometric catalog [NGAL] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   maggies - output maggies [5,NGAL]
;   ivarmaggies - corresponding inverse variance array [5,NGAL]  
;
; COMMENTS:
;   A minimum error of 0.02 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2015 Aug 8, Siena
;
; Copyright (C) 2015, John Moustakas
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

pro cfhtls_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

    nobj = n_elements(cat)
    if (nobj eq 0L) then begin
       doc_library, 'cfhtls_to_maggies'
       return
    endif

    filterlist = cfhtls_filterlist()

    names = ['u','g','r','i','z']
    errnames = names+'err'
    nband = n_elements(names)
    minerrors = replicate(0.02,nband)

    weff = k_lambda_eff(filterlist=filterlist)
    kl = k_lambda(weff,/odonnell,/silent)
    glactc, cat.ra, cat.dec, 2000.0, gl, gb, 1, /deg
    ebv = dust_getval(gl,gb,/interp,/noloop)
    
    maggies = fltarr(nband,nobj)
    ivarmaggies = fltarr(nband,nobj)
    for iband = 0L, nband-1 do begin
       ftag = tag_indx(cat[0],names[iband]+'')
       utag = tag_indx(cat[0],errnames[iband])
       good = where((cat.(ftag) gt 0.0) and (cat.(ftag) lt 30.0) and $ ; note <30 mag!
         (cat.(utag) gt 0.0) and (cat.(utag) lt 5.0),ngood) ; note 5 mag!
       if (ngood ne 0) then begin
          mag = cat[good].(ftag) - kl[iband]*ebv[good]
          magerr = cat[good].(utag)
          maggies[iband,good] = 10.0^(-0.4*mag)
          notzero = where((maggies[iband,good] gt 0.0),nnotzero)
          if (nnotzero ne 0L) then ivarmaggies[iband,good[notzero]] = $
            1.0/(0.4*alog(10.0)*(maggies[iband,good[notzero]]*magerr[notzero]))^2
       endif 
    endfor 

;   k_minerror, maggies, ivarmaggies, minerrors
    
return
end
