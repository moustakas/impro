;+
; NAME:
;   IM_DMOD2D()
;
; PURPOSE:
;   Convert from distance modulus to distance in Mpc.
;
; INPUTS: 
;   dmod - input distance modulus
;
; OPTIONAL INPUTS: 
;   sigdmod - uncertainty on DMOD
;   a_lambda - foreground Galactic extinction
;
; OUTPUTS: 
;   distance - distance in Mpc
;
; OPTIONAL OUTPUTS:
;   sigdistance - uncertainty on DISTANCE in Mpc
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 May 01, UofA
;
; Copyright (C) 2003, John Moustakas
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

function im_dmod2d, dmod, sigdmod=sigdmod, sigdistance=sigdistance, A_lambda=A_lambda
    
    ndmod = n_elements(dmod)
    if ndmod eq 0L then begin
       print, 'Syntax - distance = im_dmod2d(dmod,[sigdmod=,sigdistance=,A_lambda=])'
       return, -1L
    endif
    
    if n_elements(A_lambda) eq 0L then A_lambda = 0.0 else begin ; extinction
       if n_elements(A_lambda) gt 1L and n_elements(A_lambda) ne ndmod then begin
          print, 'A_LAMBDA and DMOD must have the same number of elements.'
          return, -1L
       endif
    endelse 

    if n_elements(sigdmod) eq 0L then sigdmod = 0.0 else begin  ; error in the distance modulus
       if n_elements(sigdmod) gt 1L and n_elements(sigdmod) ne ndmod then begin
          print, 'SIGDMOD and DMOD must have the same number of elements.'
          return, -1L
       endif
    endelse 

    dmod = dmod + A_lambda
    
    distance = 10.0^((dmod-25.0)/5.0       )     ; [Mpc]
    sigdistance = alog(10)*distance*sigdmod/5.0  ; error [Mpc]

return, distance
end
