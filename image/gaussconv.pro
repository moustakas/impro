;+
; NAME:
;	GAUSSCONV
;
; PURPOSE:
;	Convolve an image with a Gaussian kernal using Fourier
;	transforms. 
;
; INPUTS:
;	img - two-dimensional image
;
; OPTIONAL INPUTS:
;	fwhm - FWHM (= 2.35*sigma) of the Gaussian kernal (default: 3
;              pixels) 
;	npix - size of the Gaussian kernal in one dimension (default:
;              5*FWHM) 
;
; OUTPUTS:
;	An image of the same dimensions as img.
;
; PROCEDURES USED:
;	CONVOLVE() (Goddard library - uses Fourier convolution)
;
; MODIFICATION HISTORY:
;	lam98nov01ucb
;	jm01may25uofa - documented and changed to Fourier-transform
;                       convolution 
;-

;function gconvolve, x, sigma, _extra=extra
;; gaussian convolution
;; taken from VDISPFIT
;
;; special case for no smoothing
;
;    if (sigma EQ 0) then return, x
;
;    ksize = round(4*sigma+1) * 2
;    xx = findgen(ksize) - ksize/2
;
;    kernel = exp(-xx^2 / (2*sigma^2))
;    kernel = kernel / total(kernel)
;
;return, convol(x,kernel,_extra=extra)
;end

function gaussconv, img, fwhm, npix

    if n_params() le 0L then message, 'Syntax - gaussconv, img, fwhm, npix'
    
    if n_elements(fwhm) eq 0 then fwhm = 3.       ; pixels
    if n_elements(npix) eq 0 then npix = 5*fwhm   ; size of kernel
    gaux = shift(dist(npix,npix),npix/2.,npix/2.)
    gaukern = exp(-0.5*gaux^2/fwhm^2)
    gaukern = gaukern/total(gaukern)              ; normalize the area

return, convolve(img,gaukern)
end
