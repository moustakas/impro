;+
; NAME:
;   IM_STRING_STATS()
;
; PURPOSE:
;   Compute statistics on an array and build a string that looks nice
;   as a legend label. 
;
; INPUTS: 
;   arr - input array of any dimension
;
; OPTIONAL INPUTS: 
;   type - type of string array to build:
;     1 - median (mean+/-sigma)
;     2 - median+/-sigma
;     3 - like TYPE=1 but with rejection
;   ndecimal - number of decimal points to allow (default 2)
;   extra - extra keywords for IM_STATS()
;
; OUTPUTS: 
;   str - string array for use with, e.g., LEGEND
;
; COMMENTS:
;   For example, in most cases ARR will be residuals of some type. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Mar 07, NYU
;
; Copyright (C) 2009, John Moustakas
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

function im_string_stats, arr, type=type, ndecimal=ndecimal, _extra=extra

    if (n_elements(arr) eq 0L) then return, ''
    if (n_elements(ndecimal) eq 0L) then ndecimal = 2
    if (n_elements(type) eq 0L) then type = 1

    nd = string(ndecimal,format='(I0)')
    ss = im_stats(arr,_extra=extra)
    case type of
       1: str = strtrim(string(ss.median,format='(F12.'+nd+')'),2)+$
         ' ('+strtrim(string(ss.mean,format='(F12.'+nd+')'),2)+'\pm'+$
         strtrim(string(ss.sigma,format='(F12.'+nd+')'),2)+')'
       2: str = strtrim(string(ss.median,format='(F12.'+nd+')'),2)+$
         '\pm'+strtrim(string(ss.sigma,format='(F12.'+nd+')'),2)
       3: str = strtrim(string(ss.median_rej,format='(F12.'+nd+')'),2)+$
         ' ('+strtrim(string(ss.mean_rej,format='(F12.'+nd+')'),2)+'\pm'+$
         strtrim(string(ss.sigma_rej,format='(F12.'+nd+')'),2)+')'
       else: message, 'STRING type not supported'
    endcase
    
return, str
end
