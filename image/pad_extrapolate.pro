;+
; NAME:
;       PAD_EXTRAPOLATE()
;
; PURPOSE:
;       Pad and extrapolate a two-dimensional spectrum in the spatial
;       dimension using a low-order Legendre polynomial fit.
;
; CALLING SEQUENCE:
;
; INPUTS:
;       image  - two dimensional input image
;       npad   - number of spatial rows to pad
;
; OPTIONAL INPUTS:
;       ncoeff - number of Legendre coefficients to fit; default is 3
;                (order 2) 
;       ny     - number of rows to fit for the extrapolation
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       padim  - padded and extrapolated image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine is useful for preparing an image for a distortion
;       correction or registration.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 Feb 20, U of A
;-

function pad_extrapolate, image, npad, ncoeff=ncoeff, ny=ny, debug=debug

    if not keyword_set(ncoeff) then ncoeff = 3 ; quadratic
    if not keyword_set(ny) then ny = 30L
    
    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    rowaxis = findgen(nrows)

    if npad eq 0L then return, image
    
    padindx = [lindgen(npad),lindgen(npad)+npad+nrows] ; indices of padded rows
    padim = fltarr(ncols,nrows+2*npad)

    padim[*,npad:nrows+npad-1L] = image ; original image

; linear extrapolation
    
;   for i = 0L, ncols-1L do padim[i,rowindx] = $
;     interpol(padim[i,npad:nrows+npad-1L],rowaxis+npad,rowindx)

; functional extrapolation.  fit to the first and last NY rows
    
    fitrows = [lindgen(ny),lindgen(ny)+nrows-ny]
    for i = 0L, ncols-1L do begin

       y = reform(image[i,fitrows])
       fitcoeff = func_fit(fitrows,y,ncoeff,function_name='flegendre',yfit=yfit)
       padim[i,padindx] = polyleg(padindx,fitcoeff)

       if keyword_set(debug) then begin
          plot, fitrows, y, ps=4, xsty=3, ysty=3
          oplot, fitrows, yfit, thick=3.0
          cc = get_kbrd(1)
       endif

    endfor

; median extrapolation?    
    
return, padim
end

