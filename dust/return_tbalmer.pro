;+
; NAME:
;       RETURN_TBALMER()
;
; PURPOSE:
;       Return the Balmer decrement corresponding to a particular
;       electron temperature. 
;
; CALLING SEQUENCE:
;       balmer = return_tbalmer(temperature,/HaHb,/HaHg,$
;          /HbHg,/HaHd,/HdHb,/HgHb)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       temperature - electron temperature [K] (default 10000) 
;
; KEYWORD PARAMETERS:
;       HaHb - return then H-alpha/H-beta decrement (default) 
;       HaHg - return then H-alpha/H-gamma decrement
;       HbHg - return then H-beta/H-gamma decrement
;       HaHd - return then H-alpha/H-delta decrement
;       HdHb - return then H-delta/H-beta decrement
;       HgHb - return then H-gamma/H-beta decrement
;
; OUTPUTS:
;       tbalmer - the Balmer decrement at TEMPERATURE 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       BALMER_TEMPERATURE()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Aug 08 - written
;       jm04jan22uofa - modified to accept the output from
;                       BALMER_TEMPERATURE(); linearly interpolate the
;                       local temperature based on the theoretical
;                       line ratios 
;       jm04nov03uofa - documented and cleaned up
;
; Copyright (C) 2003-2004, John Moustakas
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

function return_tbalmer, temperature, HaHb=HaHb, HaHg=HaHg, $
  HbHg=HbHg, HaHd=HaHd, HdHb=HdHb, HgHb=HgHb, HgHd=HgHd

    if n_elements(temperature) eq 0L then temperature = 1E4
    
    data = balmer_temperature()

    tbalmer = interpol(data.hahb,data.temperature,temperature)                           ; Ha/Hb
    if keyword_set(HaHg) then tbalmer = interpol(data.hahg,data.temperature,temperature) ; Ha/Hg
    if keyword_set(HbHg) then tbalmer = interpol(data.hbhg,data.temperature,temperature) ; Hb/Hg
    if keyword_set(HaHd) then tbalmer = interpol(data.hahd,data.temperature,temperature) ; Ha/Hd
    if keyword_set(HdHb) then tbalmer = interpol(data.hdhb,data.temperature,temperature) ; Hd/Hb
    if keyword_set(HgHb) then tbalmer = interpol(data.hghb,data.temperature,temperature) ; Hg/Hb
    if keyword_set(HgHd) then tbalmer = interpol(data.hghd,data.temperature,temperature) ; Hg/Hd
    
return, tbalmer
end
