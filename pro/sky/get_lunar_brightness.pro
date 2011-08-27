;+
; NAME:
;   GET_LUNAR_BRIGHTNESS()
;
; PURPOSE:
;   Compute the surface brightness of the Moon [mag/arcsec2] as a
;   function of lunar age.
;
; INPUTS: 
;   age - lunar age (days)
;
; OUTPUTS: 
;   lunar - output data structure with the UBVRI surface brightness 
;
; COMMENTS:
;   Based on http://www.astro.utoronto.ca/~patton/astro/mags.html
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Mar 04, U of A
;
; Copyright (C) 2004, John Moustakas
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

function get_lunar_brightness, age

    age = [0,3,7,10,14]    ; lunar age [days]
    nage = n_elements(age)
    
    lunar = {$
      U: 0.0, $
      B: 0.0, $
      V: 0.0, $
      R: 0.0, $
      I: 0.0}
    lunar = replicate(lunar,nage)

    lunar.U = [22.0,21.5,19.9,18.5,17.0]
    lunar.B = [22.7,22.4,21.6,20.7,19.5]
    lunar.V = [21.8,21.7,21.4,20.7,20.0]
    lunar.R = [20.9,20.8,20.6,20.3,19.9]
    lunar.I = [19.9,19.9,19.7,19.5,19.2]

return, lunar
end
    
