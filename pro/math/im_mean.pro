;+
; NAME:
;   IM_MEAN()
; PURPOSE:
;   Simple wrapper on DJS_ITERSTAT to compute the mean value. 
; INPUTS:
;   x - input array
; OPTIONAL INPUTS:
;   sigrej - rejection threshold
;   maxiter - maximum number of iterations
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Feb 15, NYU
;-

function im_mean, x, sigrej=sigrej, maxiter=maxiter
  djs_iterstat, x, mean=mean, sigrej=sigrej, maxiter=maxiter
return, mean
end
