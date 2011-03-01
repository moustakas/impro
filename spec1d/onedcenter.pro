pro onedcenter, x, y, xcen, ycen, dx=dx, verbose=verbose
; jm01jul13uofa
; find the intensity-weighted (partial-pixel) center of a one-dimensional array

; dx is the tolerance (maximum fractional shift in the x-center for convergence)
; xcen is the outputed center of the array (in partial pixels) and
; ycen is the interpolated maximum value of the y-array at xcen
    
    on_error, 2
    if n_elements(x) ne n_elements(y) then message, 'Dimensions of x and y do not agree!'

    if not keyword_set(dx) then dx = 1.0E-6 ; pixels

    ymax = max(y,indx)
    xcenold = x[indx]
    weights = y/total(y) ; intensity weights

    iter = 0L
    repeat begin
       xcen = total((x-xcenold)*weights)+xcenold ; intensity-weighted center
;      xcen = total((x-xcenold)*y)/total(y)+xcenold
       dxshift = abs(xcenold-xcen)
       xcenold = xcen
       iter = iter + 1L
    endrep until (dxshift lt dx) or (iter gt 25L)

    ycen = interpol(y,x,xcen) ; linear interpolation
    
    if keyword_set(verbose) then message, 'The center converged in '+strn(iter)+' iterations.', /info
    
return
end
