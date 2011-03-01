;+
; NAME: ARM_CENTROID
;       
; CATEGORY: math
;
; PURPOSE: return image centroid coordinates (1st moments)
;
; CALLING SEQUENCE: result = ARM_CENTROID(im, [xcen, ycen, aperture, /square])
;
; INPUTS: im - 2D array to be centroided
;
; OPTIONAL INPUTS: 
;   xcen     - central x coordinate of desired region 
;              (default = center)
;   ycen     - central y coordinate of desired region 
;              (default = center)
;   aperture - radius (in pixels) of desired region; if SQUARE keyword
;              set, then this is the half-length of the square
;              aperture (default = distance to nearest edge)
;       
; KEYWORD: square - use square aperture rather than circular
;
; OUTPUTS: returns x & y centroid coordinates [x,y] (zero-indexed)
; 
; ROUTINES CALLED: DIST_CIRCLE
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., December 15, 2003
;-

function arm_centroid, image, xcen, ycen, aperture, SQUARE=square

    nx = N_ELEMENTS(image[*,0]) ; number of x elements in image
    ny = N_ELEMENTS(image[0,*]) ; number of y elements in image

; define defaults

    if not KEYWORD_SET(xcen) then xcen = FLOOR(nx / 2.0)
    if not KEYWORD_SET(ycen) then ycen = FLOOR(ny / 2.0)
    if not KEYWORD_SET(aperture) then aperture = MIN([x0, y0, nx-x0-1, ny-y0-1])

    z = REFORM(image, nx*ny)

; create x and y coordinate maps of image

    y = LINDGEN(nx*ny) / nx     
    x = LINDGEN(nx*ny) - LINDGEN(nx*ny) / nx * nx

; create map of distances from desired coordinates

    if KEYWORD_SET(square) then distmap = ABS(x-xcen) > ABS(y-ycen) $
    else DIST_CIRCLE, distmap, [nx,ny], xcen, ycen      

; ignore region outside of aperture

    ignore = WHERE(distmap gt aperture, count)
    if count gt 0L then z[ignore] = 0L

; compute first moments of image

    xc = TOTAL(z * x) / TOTAL(z)
    yc = total(z * y) / TOTAL(z)

    return, [xc, yc]

end
