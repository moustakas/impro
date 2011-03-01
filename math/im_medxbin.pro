;+
; NAME:
;   IM_MEDXBIN()
;
; PURPOSE:
;   Compute various statistical quantities in bins of fixed width.
;
; INPUTS:
;   x - independent variable [NDATA]
;   y - dependent variable [NDATA]
;   binsize - scalar bin size in units of X
;
; OPTIONAL INPUTS:
;   weight - statistical weight for each point (default 1) [NDATA] 
;   minx - minimum value of X to consider [default: min(x)]
;   maxx - maximum value of X to consider [default: max(x)]  
;   minpts - minimum number of points in each bin (default 3) 
;
; KEYWORD PARAMETERS:
;   verbose - print some of the statistics to STDOUT
;
; OUTPUTS:
;   result - statistical quantities, in each bin [NBINS]
;     xbin -  bin center
;     npts -  number of points
;     medx -  median of x
;     medy -  median of y
;     meanx -  mean of x
;     meany -  mean of y
;     sigy -  standard deviation of y
;     sigymean-  error in MEANY
;     quant05 -   5% quantile
;     quant25 -  25% quantile
;     quant75 -  75% quantile
;     quant95 -  95% quantile
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only bins containing more than MINPTS are returned (and therefore
;   the x-binning may be some factor of BINSIZE for some x).
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Jan 01, U of A - roughly based on
;      C. Tremonti's code of the similar name
;   jm05jun14uofa - added VERBOSE keyword
;   jm07sep13nyu - use Blanton's weighted_quartile() and
;     allow optional weights
;   jm08mar05nyu - cleaned up the documentation
;   jm10aug06ucsd - medium-level rewrite to be smarter and faster 
;
; Copyright (C) 2005, 2007-2008, 2010, John Moustakas
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

function im_medxbin, x, y, binsize, weight=weight1, minx=minx1, $
  maxx=maxx1, minpts=minpts1, verbose=verbose

    ndata = n_elements(x)
    ny = n_elements(y)
    nweight = n_elements(weight1)

    if (ndata eq 0L) or (ny eq 0L) or (n_elements(binsize) ne 1) then begin
       doc_library, 'im_medxbin'
       return, -1
    endif
    if (ny ne ndata) then begin    
       splog, 'X and Y must have the same number of elements'
       return, -1
    endif
    if (nweight eq 0L) then weight = y*0.0+1.0 else begin
       if (ny ne nweight) then begin
          splog, 'X,Y and WEIGHT must have the same number of elements'
          return, -1
       endif
       weight = weight1
    endelse
    
    if (n_elements(minx1) eq 0) then minx = min(x) else minx = minx1;>min(x)
    if (n_elements(maxx1) eq 0) then maxx = max(x) else maxx = maxx1;<max(x)
    if (n_elements(minpts1) eq 0) then minpts = 3 else minpts = minpts1

    inf = where((finite(x) eq 0) or (finite(y) eq 0),ninf)
    if (ninf ne 0L) then begin
       splog, 'Infinities found in X and/or Y arrays!'
       return, -1
    endif
    
; number of bins    
    nbins = abs(round((maxx-minx)/float(binsize))+1L)

; initialize the working data structure    
    res_template = {$
      xbin:     -1.0,$ ; bin center
      npts:       0L,$ ; number of points per bin
      medx:     -1.0,$ ; running X median
      meanx:    -1.0,$ ; running X mean
      medy:     -1.0,$ ; running Y median
      meany:    -1.0,$ ; running Y mean
      sigy:     -1.0,$ ; standard deviation
      sigymean: -1.0,$ ; error in MEANY
      quant05:  -1.0,$ ; 5% quantile
      quant25:  -1.0,$ ; 25% quantile
      quant75:  -1.0,$ ; 75% quantile
      quant95:  -1.0}   ; 95% quantile

    qquant = [0.5,0.05,0.25,0.75,0.95]
    for ibin = 0L, nbins-1 do begin
       res1 = res_template
       xbin = ibin * binsize + minx + (binsize/2.0)
       inbin = where((x gt xbin - (binsize/2.0)) and $
         (x le xbin + (binsize/2.0)),ninbin)

       if (ninbin lt minpts) then continue
       res1.xbin = xbin
       res1.npts = ninbin

; let IM_WEIGHTED_MEAN() handle the weighted and unweighted cases       
       if (nweight eq 0L) then begin
          res1.meanx = im_weighted_mean(x[inbin],quant=0.5,wquant=medx)
          res1.meany = im_weighted_mean(y[inbin],wsigma=sigy,$
            wmean_err=sigymean,wquant=quanty,quant=qquant)
       endif else begin
          res1.meanx = im_weighted_mean(x[inbin],$
            weights=weight[inbin],quant=0.5,wquant=medx)
          res1.meany = im_weighted_mean(y[inbin],weights=weight[inbin],$
            wsigma=sigy,wmean_err=sigymean,wquant=quanty,quant=qquant)
       endelse
       res1.medx = medx
       res1.medy = quanty[0]
       res1.sigy = sigy
       res1.sigymean = sigymean
       res1.quant05 = quanty[1]
       res1.quant25 = quanty[2]
       res1.quant75 = quanty[3]
       res1.quant95 = quanty[4]
       if (n_elements(res) eq 0L) then res = res1 else res = [temporary(res),res1]
    endfor

    if (n_elements(res) eq 0L) then begin
       splog, 'No statistics measured!  Check your MINX and MAXX inputs.'
       return, -1
    endif
    if keyword_set(verbose) then struct_print, res
    
return, res
end

