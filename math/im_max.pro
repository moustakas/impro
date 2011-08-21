;+
; NAME:
;   IM_MAX()
; PURPOSE:
;   Simple wrapper on DJS_ITERSTAT to compute a maximum value. 
; INPUTS:
;   x1 - input array
; OPTIONAL INPUTS:
;   sigrej - rejection threshold
;   maxiter - maximum number of iterations
; KEYWORD PARAMETERS:
;   nan - ignore NAN's
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Jul 13, NYU
;   jm08aug01nyu - added NAN keyword
;-

function im_max, x1, sigrej=sigrej, maxiter=maxiter, nan=nan

    if keyword_set(nan) then x = x1[where(finite(x1))] else x = x1
    djs_iterstat, x, mask=mask, sigrej=sigrej, maxiter=maxiter
    good = where(mask,ngood)
    if (ngood eq 0L) then begin
       splog, 'All points rejected!'
       return, max(x1)
    endif

return, max(x[good])
end
