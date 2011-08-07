;+
; NAME:
;   IM_EVOLVE()
; PURPOSE:
;   Compute a simple evolution model to apply to a measured absolute
;   magnitude. 
; INPUTS: 
;   z - redshift
;   q0, qz0 - parameters of the model of the form:
;     absm_evolved = absm + q0*(z-qz0)
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Feb 24, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function im_evolve, z, q0, qz0
    if (n_params() ne 3) then begin
       doc_library, 'im_evolve'
       return, -1
    endif
return, q0*(z-qz0)
end
