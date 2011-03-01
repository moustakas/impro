function gconvolve, x, sigma, _extra=extra
; gaussian convolution
; taken from VDISPFIT

; special case for no smoothing

    if (sigma EQ 0) then return, x

    ksize = round(4*sigma+1) * 2
    xx = findgen(ksize) - ksize/2

    kernel = exp(-xx^2 / (2*sigma^2))
    kernel = kernel / total(kernel)

return, convol(x,kernel,_extra=extra)
end
