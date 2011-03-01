;+
; NAME:
;   IM_LF_VMAX()
;
; PURPOSE:
;   Compute a Vmax-weighted luminosity function.
;
; INPUTS: 
;   absmag - absolute magnitude for each galaxy [NGAL]
;   oneovervmax - 1/Vmax for each object (h^3 Mpc^-3) [NGAL] 
;
; OPTIONAL INPUTS: 
;   maxabsmag - maximum absolute magnitude; see FULLBIN (default
;     max(absmag))  
;   binsize - logarithmic width of each absolute magnitude bin 
;     (default 0.1 dex)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;   binabsmag - stellar absmag bin [NBINS]
;   phi - number density (h^3 Mpc^-3) [NBINS] 
;   errphi - error on PHI (h^3 Mpc^-3) [NBINS]
;   fullbin - Boolean array indicating whether all the galaxies in
;     each bin have ABSMAG>MAXABSMAG [NBINS]
;   number - number of objects in each bin [NBINS]
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 18, UCSD 
;
; Copyright (C) 2010, John Moustakas
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

function im_lf_vmax, absmag, oneovervmax, maxabsmag=maxabsmag, $
  binsize=binsize, histmin=histmin, histmax=histmax, debug=debug

    ngal = n_elements(absmag)
    if (ngal eq 0L) then begin
       doc_library, 'im_lf_vmax'
       return, -1
    endif
    if (ngal ne n_elements(oneovervmax)) then begin
       splog, 'Dimensions of ABSMAG and ONEOVERVMAX do not agree'
       return, -1
    endif
    zero = where(oneovervmax le 0.0,nzero)
    if (nzero ne 0L) then message, 'Some ONEOVERVMAX are <= 0!'

; round the bin centers to two decimal points
    if (n_elements(binsize) eq 0) then binsize = 0.1 ; mag
    if (n_elements(maxabsmag) eq 0) then maxabsmag = max(absmag)
    if (n_elements(histmax) eq 0) then $
      histmax = ceil((max(absmag)+binsize/2.0)*100.0)/100.0-binsize/2.0

; weighted histogram    
    phi = im_hist1d(absmag,oneovervmax,binsize=binsize,obin=binabsmag,$
      binedge=0,omin=omin,omax=omax,h_err=phierr,histmin=histmin,$
      histmax=histmax,reverse_indices=rev)

; compute FULLBIN and NUMBER
    nbins = n_elements(phi)
    fullbin = intarr(nbins)+1 ; start with all full
    number = intarr(nbins)
    for ii = 0L, nbins-1 do begin
       number[ii] = rev[ii+1]-rev[ii] ; number of objects per bin
       if (rev[ii] eq rev[ii+1]) then fullbin[ii] = 0 else $
         fullbin[ii] = total(absmag[rev[rev[ii]:rev[ii+1]-1]] gt maxabsmag) eq 0.0
    endfor

; pack into a structure
    null = fltarr(50>n_elements(phi))-999.0
    lf_data = {$
      ngal:         ngal,$
      nbins:       nbins,$
      fullbin: fix(null),$
      number:  fix(null),$
      absmag:       null,$
      phi:          null,$
      phierr:       null}
    lf_data.fullbin[0:nbins-1] = fullbin
    lf_data.number[0:nbins-1] = number
    lf_data.absmag[0:nbins-1] = binabsmag
    lf_data.phi[0:nbins-1] = phi
    lf_data.phierr[0:nbins-1] = phierr

; QAplot    
    if keyword_set(debug) then begin
       ploterror, binabsmag, alog10(phi), phierr/phi/alog(10.0), $
         ps=10, xsty=3, ysty=3
       niceprint, binabsmag, phi
       cc = get_kbrd(1)
    endif

return, lf_data
end
