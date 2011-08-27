;+
; NAME:
;   COVAR2ELLIPSE()
;
; PURPOSE:
;   Given a 2D covariance matrix, compute the parameters of the
;   confidence ellipse.
;
; INPUTS: 
;   covar - 2D covariance matrix
;
; OPTIONAL INPUTS: 
;   nsigma - desired confidence interval; default is 1-sigma 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   ellipse - data structure containing the parameters of the ellipse 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Based entirely on D. Coe's beautiful Fisher-matrix
;   write-up, astro-ph/0906.3123. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jul 09, UCSD
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

function covar2ellipse, covar, nsigma=nsigma

    if (n_elements(nsigma) eq 0) then nsigma = 1.0
    
    varx = covar[0,0]
    vary = covar[1,1]
    varxy = covar[0,1]*covar[1,0]

; (squared) semi-major and semi-minor axis lengths
    aa = (varx + vary)/2.0 + sqrt(((varx - vary)^2)/4.0 + varxy)
    bb = (varx + vary)/2.0 - sqrt(((varx - vary)^2)/4.0 + varxy)

; note: the angle, theta, is defined to be positive counter-clockwise
; with respect to the x-axis, such that the parameters of the ellipse
; can be passed directly to TVELLIPSE
    theta = 90-!radeg*atan(2*sqrt(varxy)/(varx-vary))/2.0
    
; the axis lengths have to be multiplied by a coefficient that depends
; on the desired confidence level (see Table 1 in Coe's paper)
    factor = sqrt(mpchilim(errorf(nsigma/sqrt(2)),2))
    ell = {major: factor*sqrt(aa), minor: factor*sqrt(bb), angle: theta}
    
return, ell
end
