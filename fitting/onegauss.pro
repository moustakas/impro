function onegauss, xval, pp, sigmares=sigmares
; the sigma line-width is comprised of the fixed spectral resolution
; width and the variable intrinsic line width

    sigma_squared = sigmares^2.0 + pp[2]^2.0 ; total Gaussian sigma width [log-Angstrom]
    sigma = sqrt(sigma_squared)

    term1 = exp(-(xval - pp[1])^2 / (2. * sigma_squared ))
    yval = pp[0] * term1 / (sqrt(2.*!pi) * sigma)

;   print, sigmares*3E5*alog(10), pp[2]*3E5*alog(10), sigma*3E5*alog(10)
    
;   term1 = exp( - (xval - pp[1])^2 / (2. * pp[2]^2) )
;   yval = pp[0] * term1 / (sqrt(2.*!pi) * pp[2])

return, yval
end
