;+
; NAME: arm_getindex
;       
; CATEGORY: miscellaneous
;
; PURPOSE: return array index for desired element(s)
;
; CALLING SEQUENCE:
;
; INPUTS:
;   array - 1 or 2 dimensional array of numerical values
;   value - desired value, can be a single value or array
;       
; OPTIONAL INPUTS:
;
; KEYWORDS:
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;   xindex - x indices in case of 2D input array
;   yindex - y indices in case of 2D input array
;   nmatch - number of matches for each VALUE entry
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; COMMENTS: 1) double precision recommended
;           2) if VALUE is an array, then the first matching index is
;           returned for each VALUE entry; otherwise, all matching
;           indices are returned
; 
; BUG REPORT: Please report any bugs to Andrew R. Marble.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs.
;    modified to return all matching indices, ARM, 2002 June
;-

function arm_getindex, array, value, xindex=xindex, yindex=yindex, nmatch=nmatch

    nx = N_ELEMENTS(array[*,0])
    ny = N_ELEMENTS(array[0,*])
    nv = N_ELEMENTS(value)

    index  = LONARR(nv)
    xindex = LONARR(nv)
    yindex = LONARR(nv)
    nmatch = LONARR(nv)

    for i = 0, nv-1 do begin

       diff = ABS(array-value[i])
       ndx = WHERE(diff eq MIN(diff), count)
       nmatch[i] = count

       if nv gt 1 then begin

          index[i]  = ndx[0]
          yindex[i] = ndx[0] / FIX(nx)
          xindex[i] = ndx[0] - yindex[i] * nx

       endif else begin

          index  = ndx
          yindex = ndx / FIX(nx)
          xindex = ndx - yindex * nx

       endelse 

    endfor
    
;;  if nx eq 1 or ny eq 1 then begin
;;     
;;     xindex = -1L
;;     yindex = -1L
;;     
;;  endif
    
    if N_ELEMENTS(index) eq 1L then begin
       
       index  = index[0]
       xindex = xindex[0]
       yindex = yindex[0]
       
    endif 
    
    return, index

 end

