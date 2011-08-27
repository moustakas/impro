;+
; NAME:
;	GRAPHSELECT
; PURPOSE:
;	Select points within an area of a graph defined by the cursor.
;
; CALLING SEQUENCE:
;	GRAPHSELECT, x, y, indx, _extra=extra
;
; INPUTS:
;	x: the array of x-values on the plot.
;	y: the array of y-values on the plot.
;
; OPTIONAL INPUTS:
;       extra - keywords for DEFROI
;
; OUTPUTS:
;       indx - the indices of x and y that lie within the selected
;              area.  indx is equivalent to what you get with the
;              "where" function.
;
; RESTRICTIONS:
;       You cannot use this on any plot except for the most recently
;       defined plot. That is, if you made a plot (nr 1) and then
;       another one (nr 2), you can use it on nr 2 but not on nr 1. 
;
; EXAMPLE:
;	Define x=findgen(100) and y=2*x; plot, x, y. Type 
;
;	GRAPHSELECT, X, Y, INDX, /nofill
;
;	and this routine prompts you to define an area with the cursor.
;	When you are finished, the valuies of INDX are returned. You can
;	overplot the selected points to check it:
;
;	OPLOT, x[indx], y[indx], psym=2, color=128
;
; RELATED PROCEDURES:
;	IDL'S DEFROI
;
; MODIFICATION HISTORY:
;	Written by Carl Heiles. 12 Sep 1998.
;	jm02may20uofa - error checking
;       jm03jul24uofa - added EXTRA keyword
;-

pro graphselect, x, y, indx, _extra=extra

    nx = n_elements(x)
    if nx eq 0L then begin
       print, 'Syntax - graphselect, x, y, indx, _extra=extra'
       return
    endif

    if nx ne n_elements(y) then message, 'Dimensions of X and Y do not agree.'
    if !d.window[0] eq -1L then message, 'No currently defined window.'
    
; FIRST, USE DEFROI TO DEFINE THE POINTS OF INTEREST ON THE PLOT...

    result = defroi(!d.x_size,!d.y_size,_extra=extra)

; MAKE A FAKE IMAGE THAT EQUALS ZERO EVERYWHERE BUT IN THE SELECTED
; AREA...  TESTIMG IS NONZERO IN THE AREA; THE AREA IS DEFINED IN
; **DEVICE** COORDS.

    testimg = bytarr(!d.x_size,!d.y_size)
    testimg[result] = 1b

; NOW CONVERT THE DATA POINTS, WHICH ARE OF COURSE IN **DATA**
; COORDINATES, ; INTO DATA POINTS DEFINED IN **DEVICE** COORDINATES 

    xyconv = convert_coord(x,y,/data,/to_device)
    xconv = xyconv[0,*]
    yconv = xyconv[1,*]

; THEN SELECT THOSE POINTS THAT LIE WITHIN THE SELECTED REGION.

    indx = where(testimg[xconv,yconv] ne 0b)

return
end
