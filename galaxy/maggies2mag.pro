;+
; NAME:
;   MAG2MAGGIES()
;
; PURPOSE:
;   Convert from maggies to AB magnitude.
;
; INPUTS: 
;   maggies - input maggies [NBAND,NOBJ]
;
; OPTIONAL INPUTS: 
;   ivarmaggies - corresponding inverse variance maggies [NBAND,NOBJ]  
;
; OUTPUTS: 
;   mag - output AB magnitudes [NBAND,NOBJ]
;
; OPTIONAL OUTPUTS:
;   magerr - corresponding 1-sigma error [NBAND,NOBJ]
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 May 24, NYU
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

function maggies2mag, maggies, ivarmaggies=ivarmaggies, magerr=magerr
    mag = maggies*0.0
    good = where(maggies gt 0.0,ngood)
    if (ngood ne 0L) then mag[good] = -2.5*alog10(maggies[good])
    if arg_present(magerr) and (n_elements(ivarmaggies) ne 0) then begin
       magerr = ivarmaggies*0.0
       good = where((maggies gt 0.0) and (ivarmaggies gt 0.0),ngood)
       if (ngood ne 0L) then magerr[good] = 1.0/(0.4*sqrt(ivarmaggies[good])*$
         maggies[good]*alog(10.0))
    endif
return, mag
end
    
