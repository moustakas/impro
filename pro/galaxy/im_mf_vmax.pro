;+
; NAME:
;   IM_MF_VMAX()
;
; PURPOSE:
;   Compute a 1/Vmax-weighted mass *or* luminosity function.
;
; INPUTS: 
;   mass - stellar mass or absolute magnitude for each galaxy [NGAL]
;   oneovervmax - 1/Vmax for each object (h^3 Mpc^-3) [NGAL] 
;
; OPTIONAL INPUTS: 
;   binsize - logarithmic width of each stellar mass bin 
;     (default 0.1 dex)
;   masslimit - stellar mass above which the sample is complete;
;     affects LIMIT output (default min(mass)); can also be an NGAL
;     vector, which is useful when combining fields
;   minmass - minimum stellar mass
;   maxmass - maximum stellar mass
;   noutbins - number of bins in the output data structure; allows
;     different output structures from this routine to be 'stacked'  
;
; KEYWORD PARAMETERS: 
;   lf - construct the 1/Vmax-weighted luminosity function 
;   debug - make a simple plot and wait for a keystroke
;
; OUTPUTS: 
;   mf_data - the results
;     ngal - total number of galaxies
;     nbins - number of mass bins
;     number - number of objects in each bin [NBINS]
;     neff - effective number of objects in each bin [NBINS]
;     neff_err - approximate 1-sigma error in NEFF [NBINS]
;     weff - effective weight (see Zhu+09, eq 3)
;     mass - stellar mass at the center of each mass bin [NBINS]
;     phi - number density (h^3 Mpc^-3) [NBINS] 
;     phierr_poisson - simple Poisson error on PHI (h^3 Mpc^-3) [NBINS]
;     phierr_lower - lower 1-sigma error on PHI (h^3 Mpc^-3) [NBINS]
;     phierr_upper - upper 1-sigma error on PHI (h^3 Mpc^-3) [NBINS]
;     phierr - average of PHIERR_LOWER and PHIERR_UPPER (h^3 Mpc^-3) [NBINS]
;     phierr_cv - uncertainty due to cosmic variance (not computed!)
;     limit - 0=lower limit or no galaxies in the bin; 1=all points
;       have MASS>MASSLIMIT 
;
; COMMENTS:
;   Basically a glorified wrapper on IM_HIST1D().
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Aug 25, UCSD
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

function im_mf_vmax, mass, oneovervmax, binsize=binsize, masslimit=masslimit1, $
  minmass=minmass, maxmass=maxmass, noutbins=noutbins, rev=rev, lf=lf, debug=debug

    ngal = n_elements(mass)
    if (ngal eq 0L) then begin
       doc_library, 'im_mf_vmax'
       return, -1
    endif
    if (ngal ne n_elements(oneovervmax)) then begin
       splog, 'Dimensions of MASS and ONEOVERVMAX do not agree'
       return, -1
    endif
    zero = where(oneovervmax le 0.0,nzero)
    if (nzero ne 0L) then message, 'Some ONEOVERVMAX are <= 0!'

; round the bin centers to two decimal points
    if (n_elements(binsize) eq 0) then binsize = 0.1 ; dex
    nmasslimit = n_elements(masslimit1)
    if (nmasslimit eq 0) then masslimit = replicate(min(mass),ngal) else begin
       if (nmasslimit eq 1) then masslimit = replicate(masslimit1,ngal) else $
         masslimit = masslimit1
       if (n_elements(masslimit) ne ngal) then message, $
         'Dimensions of MASS and MASSLIMIT do not agree!'
    endelse 
    
    if (n_elements(minmass) eq 0) then $
      minmass = ceil((min(mass)+binsize/2.0)*100.0)/100.0-binsize/2.0

; generate weighted histogram    
    phi = im_hist1d(mass,oneovervmax,binsize=binsize,obin=binmass,$
      binedge=0,omin=omin,omax=omax,h_err=phierr_poisson,$
      histmin=minmass,histmax=maxmass,reverse_indices=rev)
    phi = phi/binsize
    phierr_poisson = phierr_poisson/binsize

; compute LIMIT, NUMBER, and the statistical error in each bin, taking
; into account small-number statistics (see Zhu+09 and Gehrels+86);
; note: the Poisson error is equivalent to sqrt(Neff)*Weff
    nbins = n_elements(phi)
    limit = intarr(nbins)+1     ; start with all complete
    number = intarr(nbins)
    weff = fltarr(nbins)
    neff = fltarr(nbins)
    neff_err = fltarr(nbins)
    medmass = fltarr(nbins)
    meanmass = fltarr(nbins)
    phierr_gehrels = fltarr(2,nbins)
    for ii = 0L, nbins-1 do begin
       number[ii] = rev[ii+1]-rev[ii] ; number of objects per bin
       if (rev[ii] eq rev[ii+1]) then begin ; NUMBER=0; compute upper limit
          limit[ii] = 0 ; no galaxies
       endif else begin
          if keyword_set(lf) then $
            limit[ii] = (binmass[ii]-binsize/2.0) le min(masslimit[rev[rev[ii]:rev[ii+1]-1]]) else $
              limit[ii] = (binmass[ii]-binsize/2.0) ge min(masslimit[rev[rev[ii]:rev[ii+1]-1]])
;         limit[ii] = total(mass[rev[rev[ii]:rev[ii+1]-1]] lt masslimit) eq 0.0
          weff[ii] = total(oneovervmax[rev[rev[ii]:rev[ii+1]-1]]^2,/double)/$ ; effective weight
            total(oneovervmax[rev[rev[ii]:rev[ii+1]-1]],/double)
          if (weff[ii] le 0) then message, 'This should not happen!'
          neff[ii] = total(oneovervmax[rev[rev[ii]:rev[ii+1]-1]],/double)/weff[ii] ; effective number of objects
          neff_err1 = im_poisson_limits(neff[ii],0.8413) ; 1-sigma
          neff_err[ii] = (neff_err1[1]-neff_err1[0])/2.0 ; error is asymmetric, but approximate
          phierr_gehrels[0,ii] = phi[ii]-weff[ii]*neff_err1[0]/binsize ; lower
          phierr_gehrels[1,ii] = weff[ii]*neff_err1[1]/binsize-phi[ii] ; upper
          medmass[ii] = djs_median(mass[rev[rev[ii]:rev[ii+1]-1]]) 
          meanmass[ii] = djs_mean(mass[rev[rev[ii]:rev[ii+1]-1]]) 
       endelse
    endfor

; pack into a structure
    if (n_elements(noutbins) eq 0) then noutbins = nbins
    null = fltarr(noutbins)

    if keyword_set(lf) then begin
       mf_data = {$
         ngal:           ngal,$
         binsize:     binsize,$
         nbins:         nbins,$
         number:    fix(null),$
         neff:           null,$
         neff_err:       null,$
         weff:           null,$
         absmag:         null,$
         phi:            null,$
         phierr_poisson: null,$
         phierr_lower:   null,$
         phierr_upper:   null,$
         phierr:         null,$
         phierr_cv:      null,$
         limit:      fix(null)}
       mf_data.absmag[0:nbins-1] = binmass
    endif else begin
       mf_data = {$
         ngal:           ngal,$
         binsize:     binsize,$
         nbins:         nbins,$
         number:    fix(null),$
         neff:           null,$
         neff_err:       null,$
         weff:           null,$
         mass:           null,$
         medmass:        null,$
         meanmass:       null,$
         phi:            null,$
         phierr_poisson: null,$
         phierr_lower:   null,$
         phierr_upper:   null,$
         phierr:         null,$
         phierr_cv:      null,$
         limit:      fix(null)}
       mf_data.mass[0:nbins-1] = binmass
    endelse

    mf_data.limit[0:nbins-1] = limit
    mf_data.number[0:nbins-1] = number
    mf_data.neff[0:nbins-1] = neff
    mf_data.neff_err[0:nbins-1] = neff_err
    mf_data.weff[0:nbins-1] = weff
    mf_data.medmass[0:nbins-1] = medmass
    mf_data.meanmass[0:nbins-1] = meanmass
    
    mf_data.phi[0:nbins-1] = phi
    mf_data.phierr_poisson[0:nbins-1] = phierr_poisson
    mf_data.phierr_lower[0:nbins-1] = reform(phierr_gehrels[0,*])
    mf_data.phierr_upper[0:nbins-1] = reform(phierr_gehrels[1,*])
    mf_data.phierr[0:nbins-1] = total(phierr_gehrels,1)/2.0 ; average error
       
; QAplot    
    if keyword_set(debug) then begin
       gd = where(phi gt 0)
       phierr = phierr_poisson[gd]/phi[gd]/alog(10.0)
       ploterror, binmass[gd], alog10(phi[gd]), phierr, ps=10, xsty=3, ysty=3
       niceprint, binmass, limit, number, neff, weff, phi, phierr_poisson
    endif

return, mf_data
end
