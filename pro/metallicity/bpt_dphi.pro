;+
; NAME:
;   BPT_DPHI()
;
; PURPOSE:
;   Compute the BPT classification "D" and "phi" parameters of
;   Kauffmann+03. 
;
; INPUTS: 
;   niiha - logarithmic [NII]/Ha ratio
;   oiiihb - logarithmic [OIII]/Hb ratio
;
; OUTPUTS: 
;   d - Kauffmann's "D" parameter
;   phi - Kauffmann's "phi" parameter
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Feb 24, NYU
;
; Copyright (C) 2008, John Moustakas
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

function bpt_dphi, niiha, oiiihb, phi=phi

    if (n_elements(niiha) eq 0L) or (n_elements(oiiihb) eq 0L) then begin
       doc_library, 'bpt_dphi'
       return, -1L
    endif
    
    x = niiha + 0.45
    y = oiiihb + 0.5
    d = sqrt(x^2 + y^2)
    phi = atan(x/y)*!radeg

return, d
end
