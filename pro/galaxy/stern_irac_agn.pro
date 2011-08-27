;+
; NAME:
;   STERN_IRAC_AGN()
;
; PURPOSE:
;   Identify AGN using the Stern wedge.
;
; INPUTS: 
;   ch1_ab, ch2_ab, ch3_ab, ch4_ab - IRAC ch1-4 photometry (AB mag) 
;
; OUTPUTS: 
;   isagn - Boolean array (1=AGN)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 16, UCSD
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

function stern_irac_agn, ch1_ab, ch2_ab, ch3_ab, ch4_ab

    v2ab = k_vega2ab(filterlist=irac_filterlist(),/kurucz,/silent)
    ch12 = (ch1_ab-ch2_ab)-(v2ab[0]-v2ab[1])
    ch34 = (ch3_ab-ch4_ab)-(v2ab[2]-v2ab[3])

    isagn = (ch34 gt 0.6) and (ch12 gt (0.2*ch34+0.18)) and $
      (ch12 gt (2.5*ch34-3.5))
    
return, isagn
end
