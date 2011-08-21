;+
; NAME:
;   ONEDCENTER()
;
; PURPOSE:
;   Find the intensity-weighted (partial-pixel) center of a 1D array.  
;
; INPUTS: 
;   x - independent input vector (e.g., wavelength) [NPTS]
;   y - dependent input vector (e.g., flux) [NPTS]
;
; OPTIONAL INPUTS: 
;   dx - tolerance; maximum fractional shift in the x-center for
;     convergence 
;
; KEYWORD PARAMETERS: 
;   verbose - print messages to STDOUT
;
; OUTPUTS: 
;   xcen - partial-pixel center of X
;   ycen - interpolated value of Y at XCEN
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Jul 13, U of A
;
; Copyright (C) 2001, John Moustakas
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

pro onedcenter, x, y, xcen, ycen, dx=dx, verbose=verbose

    npts = n_elements(x)
    if (npts eq 0L) or (npts ne n_elements(y)) then begin
       doc_library, 'onedcenter'
       return
    endif

    if (n_elements(dx) eq 0) then dx = 1.0E-6 ; pixels

    ymax = max(y,indx)
    xcenold = x[indx]
    weights = y/total(y) ; intensity weights

    iter = 0L
    repeat begin
       xcen = total((x-xcenold)*weights)+xcenold ; intensity-weighted center
;      xcen = total((x-xcenold)*y)/total(y)+xcenold
       dxshift = abs(xcenold-xcen)
       xcenold = xcen
       iter++
    endrep until (dxshift lt dx) or (iter gt 25)

    ycen = interpol(y,x,xcen) ; linear interpolation
    
    if keyword_set(verbose) then message, 'The center converged in '+$
      strn(iter)+' iterations.', /info
    
return
end
