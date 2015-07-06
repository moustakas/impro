;+
; NAME:
;   UNWISE_TO_MAGGIES
;
; PURPOSE:
;   Convert the Tractor unWISE photometry to AB maggies.  
;
; INPUTS: 
;   cat - input photometric catalog [NGAL] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   maggies - output maggies [2,NGAL]
;   ivarmaggies - corresponding inverse variance array [8,NGAL]  
;
; COMMENTS:
;   Assume no minimum magnitude uncertainty.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 May 14, Siena
;
; Copyright (C) 2014, John Moustakas
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

pro unwise_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, $
  prefix=prefix, ratag=ratag, dectag=dectag, nodust=nodust, ebv=ebv

    ngal = n_elements(cat)
    if (ngal eq 0L) then begin
       doc_library, 'unwise_to_maggies'
       return
    endif

    filterlist = wise_filterlist(/short)
    weff = k_lambda_eff(filterlist=filterlist)
    nbands = n_elements(filterlist)

    if n_elements(prefix) eq 0 then prefix2 = '' else prefix2 = prefix+'_'
    tags = prefix2+['w1','w2']+'_nanomaggies'
    ivartags = prefix2+['w1','w2']+'_nanomaggies_ivar'

; http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
    v2ab = [2.699,3.339]
    kl = ext_ccm(weff)*3.1
;   kl = k_lambda(weff,/ccm,/silent)
    
; correct for dust
    if keyword_set(nodust) eq 0 then begin
       if n_elements(ebv) ne 0L then begin
          if n_elements(ebv) ne ngal then message, 'Dimensions of CAT and EBV must match!'
       endif else begin
          if n_elements(ratag) eq 0 then ratag = 'ra'
          if n_elements(dectag) eq 0 then dectag = 'dec'
          glactc, cat.(tag_indx(cat,ratag)), cat.(tag_indx(cat,ratag)), $
            2000.0, gl, gb, 1, /deg
          ebv = dust_getval(gl,gb,/interp,/noloop)
       endelse
    endif else ebv = 0.0
    
; conversion factor from nanomaggies to maggies
    fact = 1D-9

; construct maggies and ivarmaggies in each band       
    maggies = dblarr(nbands,ngal)
    ivarmaggies = dblarr(nbands,ngal)
    for ib = 0, nbands-1 do begin
       ftag = tag_indx(cat[0],tags[ib])
       itag = tag_indx(cat[0],ivartags[ib])
       maggies[ib,*] = cat[*].(ftag)*fact*10D^(0.4*(kl[ib]*ebv-v2ab[ib]))
       ivarmaggies[ib,*] = cat[*].(itag)/(fact*10D^(0.4*(kl[ib]*ebv-v2ab[ib])))^2D
    endfor

return   
end
