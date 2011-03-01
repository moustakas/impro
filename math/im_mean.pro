; wrapper for djs_iterstat
; J. Moustakas, 2008 Feb 15, NYU

function im_mean, x, sigrej=sigrej, maxiter=maxiter

  djs_iterstat, x, mean=mean, sigrej=sigrej, maxiter=maxiter

  return, mean
end
