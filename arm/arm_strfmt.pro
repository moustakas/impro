;+
; NAME: ARM_STRFMT
;       
; CATEGORY: string manipulation
;
; PURPOSE: make all string array elements the same length
;
; CALLING SEQUENCE: ARM_STRFMT, array, [padchar=, truncate=, length=, \reverse]
;
; INPUTS: array - array of string elements
;       
; OPTIONAL INPUTS:
;    padchar  - character to pad array elements with 
;               (default = blank space)
;    position - position (from end or beginning depending on if the
;               reverse keyword is set) within string to place padding
;               (default = 0)
;    truncate - character to truncate array elements with 
;               (default = notruncation)
;
; KEYWORDS: reverse - pad beginning of array elements instead of end
;
; OPTIONAL OUTPUTS: length - new length of array elements
;
; EXAMPLE: 
;    IDL> x = ARM_STRFMT(['abc','a','abcd','ab'], padchar='*')
;
;    x -> abc* a*** abcd ab**
;
;    IDL> x = ARM_STRFMT(['abc','a','abcd','ab'], padchar='*', position=1)
;
;    x -> a*bc a*** abcd a**b
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., August 2001
;    padchar and reverse features added by A.R.Marble, October 2003
;-

function ARM_STRFMT, array, padchar=padchar, truncate=truncate, $
            position=position, length=length, reverse=reverse
    
    if not KEYWORD_SET(padchar) then padchar = ' '

    n = N_ELEMENTS(array)
    l = INTARR(n)

; trim leading and trailing blank spaces

    copy = STRTRIM(array, 2)

; truncate array elements if desired

    if N_ELEMENTS(truncate) ne 0 then for i=0,n-1 do $
      copy[i] = (STRSPLIT(copy[i], truncate, /extract))[0]
      
    if N_ELEMENTS(position) ne 0 then begin

       if KEYWORD_SET(reverse) then begin
          
          excise = STRMID(copy, 0, position)
          copy = STRMID(copy, position)
          
       endif else begin
          
          excise = STRMID(copy, position-1, /reverse_offset)

          for i=0,n-1 do copy[i] = $
            STRMID(copy[i], 0, STRPOS(copy[i], excise[i], /reverse_search))
          
       endelse

    endif

; determined lengths of array elements

    for i=0,n-1 do l[i] = STRLEN(copy[i])
    
    if N_ELEMENTS(length) eq 0 then length = max(l)             ; desired string length

    for i=0,n-1 do begin

       diff = length - l[i]     ; length difference

       if diff gt 0 then begin

; pad the beginning or end of string

          for j=0,diff-1 do if KEYWORD_SET(reverse) then $
            copy[i] = padchar + copy[i] else $
            copy[i] = copy[i] + padchar

       endif

    endfor
    
    if N_ELEMENTS(position) ne 0 then begin
    
       if KEYWORD_SET(reverse) then copy = excise + copy $
       else copy = copy + excise

    endif

    return, copy
    
 end





