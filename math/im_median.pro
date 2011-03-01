; wrapper for djs_iterstat
; J. Moustakas, 2008 Feb 15, NYU

function im_median, x, sigrej=sigrej, maxiter=maxiter

  djs_iterstat, x, median=median, sigrej=sigrej, maxiter=maxiter

  return, median
end
