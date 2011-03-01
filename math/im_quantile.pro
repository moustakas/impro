;+
; NAME:
;   IM_QUANTILE()
;
; PURPOSE:
;   Given a set of values and weights, compute the weighted
;   quantile(s).  
;
; INPUTS: 
;   values - input values
;   weights - corresponding weights (assumed to be unity if not given) 
;
; OPTIONAL INPUTS: 
;   quant - desired quantiles (default 0.5, i.e., the median)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   quantile - the weighted quantile of the input array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is analogous to M. Blanton's WEIGHTED_QUANTILE() but
;   it uses interpolation rather than binning, making it more
;   accurate, but slower.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 05, UCSD
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

function im_quantile, values, weights, avg=avg, quant=quant

    message, 'This is wrong!'
    
    npts = n_elements(values)
    if (npts eq 0L) then begin
       doc_library, 'im_quantile'
       return, -1
    endif

    if (n_elements(quant) eq 0) then quant = 0.5
    if (n_elements(weights) eq 0) then weights = values*0.0+1.0
    if (n_elements(weights) ne npts) then begin
       splog, 'Dimensions of VALUES and WEIGHTS must match'
    endif

    srt = sort(values)
    svalues = values[srt]
    sweights = weights[srt]
    scum = total(sweights,/cum,/double)/total(sweights,/double)
;   scum = total(sweights,/cum,/double)
;   scum = scum/scum[npts-1L]
    quantile = interpol(svalues,scum,quant)

    if arg_present(avg) then avg = total(values*weights,/double)/$
      total(weights,/double)

;   message, 'This routine does not work as advertised'
;   djs_plot, svalues, scum, ps=-6, xsty=3, ysty=3, yr=[0,1]
    djs_plot, svalues, scum, ps=-6, xsty=3, ysty=3, xr=[-0.1,0.1]
    djs_oplot, quantile[0]*[1,1], !y.crange, color='cyan'
    
return, quantile
end
