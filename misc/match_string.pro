;+
; NAME:
;   MATCH_STRING()
;
; PURPOSE:
;   Matches strings one letter at a time.
;
; INPUTS: 
;   instring - input string (can be a vector)
;   allstrings - list of strings againsts which to match INSTRING 
;
; KEYWORD PARAMETERS: 
;   exact - match the string exactly 
;   silent - suppress messages to STDOUT
;
; OUTPUTS: 
;   matchstring - nearest match of INSTRING in ALLSTRINGS
;
; OPTIONAL OUTPUTS: 
;   findex - index of the first match
;   index - indices of all matches
;
; COMMENTS:
;   This routine is pretty dumb and could be done trivially with
;   regular expression matching with existing routines.
;
; EXAMPLES:
;   IDL> print, match_string('bob',['bobo','jimbob','bobobo'],findex=ff,index=ii)
;     bobo
;   IDL> print, ff
;     0
;   IDL> print, ii
;     0           2
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Oct 18, U of A
;   jm02jun25uofa - added exact keyword
;   jm02nov12uofa - added FINDEX keyword (first match)    
;   jm03may22uofa - fixed case if NAME is ''
;   jm06jan19uofa - check for STRMATCH wildcard characters    
;
; Copyright (C) 2001-2003, 2006, John Moustakas
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

function match_string, instring, allstrings, index=index, $
  findex=findex, exact=exact, silent=silent
    
    ninstring = n_elements(instring)
    if ninstring gt 1L then begin
       findex = lonarr(ninstring)
       for ii = 0L, ninstring-1L do begin
          matchstring1 = match_string(instring[ii],allstrings,index=index1,$
            findex=findex1,exact=exact,silent=silent)
          if ii eq 0 then begin
             index = index1 
             matchstring = matchstring1 
             findex = findex1
          endif else begin
             index = [index,index1]
             matchstring = [matchstring,matchstring1]
             findex = [findex,findex1]
          endelse
       endfor
       return, matchstring
    endif
    
    name = strlowcase(strcompress(instring,/remove))
    name = repstr(repstr(name,'[',''),']')
    len = strlen(name)
    newname = strarr(len>1)
    for i = 0L, len-1L do newname[i] = strmid(name,i,1)

    std = strlowcase(strcompress(allstrings,/remove))
    std = repstr(repstr(std,'[',''),']')

    if keyword_set(exact) then begin

       index = where(strmatch(std,name) eq 1B,nmatch)
       if (nmatch eq 0L) then begin
          if not keyword_set(silent) then message, $
            'Exact match not found for '+strupcase(name), /info
          matchstring = ''
          findex = -1L
       endif else begin
          matchstring = allstrings[index]
          findex = index[0]
       endelse

    endif else begin
    
; cycle on each letter in the INSTRING
       j = -1L
       repeat begin
          j = j+1L
;         match = where(strcmp(strjoin(newname[0:j]),allstrings,j+1L,/fold_case) eq 1B,nmatch)
          index = where(strcmp(strjoin(newname[0:j]),std,j+1L,/fold_case) eq 1B,nmatch)
       endrep until (nmatch eq 1B) or (j eq len-1L)

       if (nmatch eq 0B) and (j eq len-1L) then matchstring = '' else begin
          matchstring = (allstrings[index])[0]
          findex = index[0]
       endelse
;      message, 'No data file matches!' else matchstring = (allstrings[index])[0]

    endelse

return, matchstring
end

