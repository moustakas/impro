;+
; NAME:
;   IM_STRSPLIT()
;
; PURPOSE:
;   Simple vectorized version of STRSPLIT().
;
; INPUTS: 
;   str - input string (vector) array
;   key - scalar pattern on which to do the splitting
;
; OPTIONAL INPUTS: 
;   _extra - extra keywords for STRSPLIT()
;
; OUTPUTS: 
;   split - parsed string
;
; COMMENTS:
;   This is a pretty dumb piece of code and highly susceptible to
;   breaking.  
;
; EXAMPLES:
;   IDL> str = ['2006-Jan-10','2009-Feb-10']
;   IDL> key = '-' then: 
;   IDL> split = im_strsplit(str,key,/extract)
;   IDL> help, split  
;   SPLIT           STRING    = Array[3, 2]
;   IDL> print, split
;   2006 Jan 10
;   2009 Feb 10
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 May 21, NYU
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

function im_strsplit, str, key, _extra=extra

    nstr = n_elements(str)
    if (nstr eq 0) or (n_elements(key) eq 0) then begin
       doc_library, 'im_strsplit'
       return, -1
    endif
    split1 = strsplit(str[0],key,_extra=extra)
    split = strarr(n_elements(split1),nstr)
    for ii = 0L, nstr-1 do split[*,ii] = strsplit(str[ii],key,_extra=extra)
    
return, split
end
    
