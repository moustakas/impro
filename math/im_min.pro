;+
; NAME:
;   IM_MIN()
; PURPOSE:
;   Simple wrapper on DJS_ITERSTAT to compute the minimum value. 
; INPUTS:
;   x1 - input array
; OPTIONAL INPUTS:
;   sigrej - rejection threshold
;   maxiter - maximum number of iterations
; KEYWORD PARAMETERS:
;   ignorezero - exclude zero from the calculation
;   nan - ignore NAN's
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Jul 13, NYU
;   jm08aug01nyu - added NAN keyword
;-

function im_min, x1, sigrej=sigrej, maxiter=maxiter, $
  ignorezero=ignorezero, nan=nan

    if keyword_set(nan) then x = x1[where(finite(x1))] else x = x1
    djs_iterstat, x, mask=mask, sigrej=sigrej, maxiter=maxiter
    if keyword_set(ignorezero) then $
      good = where((mask eq 1) and (x ne 0.0),ngood) else $
        good = where((mask eq 1),ngood)
    if (ngood eq 0L) then begin
       splog, 'All points rejected!'
       return, min(x1)
    endif

return, min(x[good])
end
