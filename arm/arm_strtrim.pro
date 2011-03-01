;+
; NAME: ARM_STRTRIM
;       
; CATEGORY: string manipulation
;
; PURPOSE: trim leading/trailing characters
;
; CALLING SEQUENCE: result = ARM_STRTRIM(str, pattern, [/opposite, /leading]) 
;
; INPUTS:
;   str     - charcter string to be trimmed
;   pattern - string of one or more characters to be trimmed from
;             beginning/end of STR (or not to be trimmed if OPPOSITE
;             keyword is set)
;       
; KEYWORDS:
;   opposite - trim anything from beginning or end other than PATTERN
;   leading  - trim leading characters rather than trailing 
;
; OUTPUTS: trimmed character string is returned
;
; EXAMPLE: IDL> a = STRCOMPRESS(1.0250000, /remove_all)
;          IDL> result = ARM_STRTRIM(a, '0')
;         
;               result -> '1.025'
;
; COMMENTS: same as STRTRIM if pattern = ' '
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., December 3, 2003
;-
function arm_strtrim, str, pattern, opposite=opposite, leading=leading

    if KEYWORD_SET(leading) then reverse_offset = 0L else reverse_offset = 1L

    done = 0L

    while not done do begin

       if not KEYWORD_SET(leading) then begin
          start1 = STRLEN(pattern)-1 
          start2 = 0L
       endif else begin 
          start1 = 0L
          start2 = STRLEN(pattern)
       endelse

       lastchar = STRMID(str, start1, STRLEN(pattern), reverse_offset=reverse_offset) 

       if (KEYWORD_SET(opposite) and lastchar ne pattern) or $
         (not KEYWORD_SET(opposite) and lastchar eq pattern) then $
         str = STRMID(str, start2, STRLEN(str)-STRLEN(pattern)) else done = 1L

    endwhile
    
return, str

end   
