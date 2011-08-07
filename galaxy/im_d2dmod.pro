;+
; NAME:
;   IM_DMOD2D()
;
; PURPOSE:
;   Convert from distance in Mpc to distance modulus.
;
; INPUTS: 
;   dist - input distance in Mpc
;
; OPTIONAL INPUTS: 
;   err_dist - uncertainty on DIST [Mpc]
;
; OUTPUTS: 
;   dmod - distance modulus
;
; OPTIONAL OUTPUTS:
;   err_dmod - uncertainty on DMOD
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Jun 06, NYU
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

function im_d2dmod, dist, err_dist=err_dist, err_dmod=err_dmod
    
    ndist = n_elements(dist)
    if (ndist eq 0L) then begin
       doc_library, 'im_d2dmod'
       return, -1L
    endif
    
    dmod = 5.0*alog10(dist)+25.0
    if arg_present(err_dmod) then begin
       if (n_elements(err_dist) eq 0L) then err_dist = dist*0.0 else begin
          if (n_elements(err_dist) ne ndist) then begin
             print, 'Dimensions of ERR_DIST and DIST must agree'
             return, -1L
          endif
       endelse
       err_dmod = 5.0*err_dist/dist/alog(10.0)
    endif 
    
return, dmod
end
