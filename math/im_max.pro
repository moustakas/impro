; wrapper for djs_iterstat
; J. Moustakas, 2008 June 13, NYU
; jm08aug01nyu - added NAN keyword

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
