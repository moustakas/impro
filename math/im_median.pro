;+
; NAME:
;   IM_MEDIAN()
; PURPOSE:
;   Simple wrapper on DJS_ITERSTAT to compute the median value. 
; INPUTS:
;   x - input array
; OPTIONAL INPUTS:
;   sigrej - rejection threshold
;   maxiter - maximum number of iterations
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Feb 15, NYU
;-

function im_median, x, sigrej=sigrej, maxiter=maxiter
  djs_iterstat, x, median=median, sigrej=sigrej, maxiter=maxiter
return, median
end
