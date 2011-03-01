;+
; NAME:
;   IM_WEIGHTED_MEAN()
;
; PURPOSE:
;   Compute the weighted mean, standard deviations, and error in the
;   mean, given an array of values and statistical uncertainties.
;
; INPUTS: 
;   values - input array of values [NPTS]
;
; OPTIONAL INPUTS: 
;   errors - uncertainty for each element of VALUES [NPTS]
;   weights - statistical weight for each element of VALUES [NPTS]; if
;     both ERRORS and WEIGHTS are passed then ERRORS takes precedence 
;   quant - see Blanton's WEIGHTED_QUANTILE() (only used if
;     WQUANTILE is desired)
;
; KEYWORD PARAMETERS: 
;   nongauss - assume non-Gaussian errors; see COMMENTS, below; if
;   using /NONGAUSS then you should pass WEIGHTS, not ERRORS
;
; OUTPUTS: 
;   wmean - weighted mean
;
; OPTIONAL OUTPUTS:
;   wsigma - weighted standard deviation
;   wmean_err - weighted error in mean
;   wquantile - weighted quantiles for each value of QUANT (see 
;     WEIGHTED_QUANTILE) 
;
; COMMENTS:
;   This routine can do several different things and is loosely based
;   on E. Sheldon's WMOM routine.  The default use would be to
;   pass VALUES and the 1-sigma ERRORS; the routine would then return
;   the weighted mean and, optionally, the weighted error in the mean
;   and weighted standard deviation.  If the errors are not
;   Gaussian-distributed, however, then you should also set /NONGAUSS,
;   in which case the weighted error in the mean is computed from the
;   distribution itself (talk to E. Sheldon for details and also see
;   various wiki math pages).  
;
;   The second behavior is to pass WEIGHTS *instead* of ERRORS.  For
;   Gaussian errors WEIGHTS=1/ERRORS^2 (i.e., the inverse variance),
;   but nevertheless this optional input can be very useful.  If both
;   ERRORS *and* WEIGHTS are passed then ERRORS takes precedence.
;
;   Note that if both ERRORS *and* WEIGHTS are missing then uniform
;   weights are assumed and the unweighted mean, unweighted error in
;   the mean, and unweighted standard deviation are returned.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Mar 28, NYU - streamlined version of
;     E. Sheldon's WMOM.PRO routine
;   jm10jul11ucsd - some small updates
;   jm10aug06ucsd - substantial rewrite and expanded documentation 
;
; Copyright (C) 2008, 2010, John Moustakas
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

function im_weighted_mean, values, errors=errors, weights=weights, $
  wsigma=wsigma, wmean_err=wmean_err, wquantile=wquantile, quant=quant, $
  nongauss=nongauss
    
    npts = n_elements(values)
    if (npts eq 0L) then begin
       doc_library, 'im_weighted_mean'
       return, -1.0
    endif

    if (n_elements(errors) eq 0L) and (n_elements(weights) eq 0L) then begin
       noweights = 1
       weights = values*0.0+1.0
    endif else noweights = 0

    if (n_elements(errors) eq 0L) then begin
       if (n_elements(weights) ne npts) then begin
          splog, 'Dimensions of VALUES and WEIGHTS do not agree'
          return, -1.0
       endif
       bad = where((finite(weights) eq 0),nbad)
       if (nbad ne 0L) then begin
          splog, 'WEIGHTS contains infinities!'
          return, -1
       endif
    endif else begin
       if (n_elements(errors) ne npts) then begin
          splog, 'Dimensions of VALUES and ERRORS do not agree'
          return, -1.0
       endif
       bad = where((errors le 0.0) or (finite(errors) eq 0),nbad)
       if (nbad ne 0L) then begin
          splog, 'ERRORS contains negative numbers or infinities!'
          return, -1
       endif
       if keyword_set(nongauss) then begin
          splog, 'Non-Gaussian weights not defined!'
          return, -1
       endif
       weights = 1.0/errors^2 ; assume Gaussian errors
    endelse

; compute the weighted mean and weighted standard deviation
    wtot = total(weights,/double)
    if (wtot eq 0.0) then begin
       splog, 'All the weights are zero!'
       return, 0.0
    endif
    wmean = total(weights*values,/double)/wtot 

; standard deviation
    if arg_present(wsigma) then begin
       if noweights then wsigma = djsig(values) else $
         wsigma = sqrt(total(weights*(values-wmean)^2)/wtot)
    endif

; error in the mean
    if arg_present(wmean_err) then begin
       if noweights then wmean_err = wsigma/sqrt(npts) else begin
          if (n_elements(errors) eq 0L) then begin                      ; no errors supplied
             wmean_err = sqrt(total(weights^2*(values-wmean)^2)/wtot^2) ; weighted error in the mean
          endif else begin
; use the formula if the errors are non-gaussian, otherwise compute
; the usual standard error on the mean
             if keyword_set(nongauss) then $
               wmean_err = sqrt(total(weights^2*(values-wmean)^2)/wtot^2) else $
                 wmean_err = 1.0/sqrt(wtot) 
          endelse
       endelse
    endif

; finally, optionally compute the weighted median
    if arg_present(wquantile) then wquantile = $
      weighted_quantile(values,weights,quant=quant)
    
return, wmean
end
