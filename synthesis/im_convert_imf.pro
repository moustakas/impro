;+
; NAME:
;       IM_CONVERT_IMF()
;
; PURPOSE:
;       Return the conversion factor between various IMFs and the
;       Salpeter IMF between 0.1-100 M_sun.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;       from_chabrier      - 
;       from_kroupa01      - 
;       from_baldry03      - 
;       from_diet_salpeter - 
;
; OUTPUTS: 
;       The value that must be *added* to convert *to* Salpeter
;       (1955), unless LINEAR=1, in which case the value returned must
;       be multiplied.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       * Chabrier (2003) - See end of Section 3.2 in Gallazzi et
;                           al. (2008).  I have also checked this
;                           value by comparing the BC03 models for
;                           each IMF, obtaining 0.25 dex.
;       * Kroupa (2001)   - See Brinchmann et al. (2004).  Note that
;                           Bell et al. (2003) gives 0.3 dex!!
;       * Baldry & Glazebrook (2003) - See also Savaglio et
;                                      al. (2005). 
;       * Diet Salpeter   - See Bell et al. (2003).
;       
;       Note that a variety of conversion factors are also given in
;       Wilkins et al. 2008.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2008, Apr 05, NYU - written
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

function im_convert_imf, from_chabrier=from_chabrier, from_kroupa01=from_kroupa01, $
  from_baldry03=from_baldry03, from_diet_salpeter=from_diet_salpeter, linear=linear
    
    if keyword_set(from_chabrier) then factor = alog10(1.75)
    if keyword_set(from_kroupa01) then factor = alog10(1.5)
    if keyword_set(from_baldry03) then factor = alog10(1.8)
    if keyword_set(from_diet_salpeter) then factor = +0.15

    if keyword_set(linear) then return, 10.0^factor else return, factor

end
