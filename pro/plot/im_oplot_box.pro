;+
; NAME:
;   IM_OPLOT_BOX
;
; PURPOSE:
;   Overplot an arbitrarily positioned (and rotated) box on an
;   existing plot.  
;
; INPUTS:
;   xsize - length of the box along the X axis
;   ysize - length of the box along the Y axis 
;   pa    - rotation angle measured counter-clockwise from the y
;           axis  (identical to astronomical position angles)
;           [degrees] 
;
; OPTIONAL INPUTS:
;   xoffset - see IM_OFFSET_AND_ROTATE()
;   yoffset - see IM_OFFSET_AND_ROTATE()
;   extra   - keywords for DJS_OPLOT
;
; KEYWORD PARAMETERS:
;   noplot  - do not make the plot
;
; OUTPUTS:
;   A box with the specified properties is plotted.
;
; OPTIONAL OUTPUTS:
;   corners - data coordinates of the box corners
;
; PROCEDURES USED:
;   IM_OFFSET_AND_ROTATE(), DJS_OPLOT
;
; COMMENTS:
;   Y is "along" the slit, X is perpendicular to the slit. 
;
; EXAMPLE:
;   Try drawing a box around an existing plot of a galaxy.  The
;   major and minor axes are 60" and 30" respectively.  An image
;   of the galaxy with coordinate axes in arcsec has already been
;   displayed.  The position angle of the galaxy is 30 degrees.
;
;   IDL> pa = 30.0
;   IDL> im_oplot_box, 30.0, 60.0, 30.0, line=0, thick=3.0
;
;   To display a simple rotated box try.
;
;   IDL> plot, [0], [0], /nodata, xrange=[-2,2], yrange=[-2,2], xsty=3, ysty=3
;   IDL> im_oplot_box, 1.0, 3.0, 30, line=0, thick=3.0
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 May 22, U of A
;   jm03july10uofa - added NOPLOT keyword
;
; Copyright (C) 2003, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro im_oplot_box, xsize, ysize, pa, xoffset=xoffset, yoffset=yoffset, $
  corners=corners, noplot=noplot, _extra=extra

    if n_params() ne 3L then begin
       print, 'Syntax - im_oplot_box, xsize, ysize, pa, [xoffset=, yoffset=],$'
       print, '   /noplot, [_extra=extra]'
       return
    endif

; establish the coordinates (corners) of the box clock-wise from the
; lower left corner:

; xy = [ [min(x),min(y)], [min(x),max(y)], [max(x),max(y)], [max(x),min(y)] ]

    xy = [ [-1,-1], [-1,1], [1,1], [1,-1] ] * ( (0.5*[xsize,ysize]) # (fltarr(4)+1))
    
    corners = im_offset_and_rotate(xy,pa,xoffset=xoffset,yoffset=yoffset)
    indx = [0,1,2,3,0]

    if not keyword_set(noplot) then begin
       djs_oplot, [corners[0,indx],corners[0,indx+1]], $
         [corners[1,indx],corners[1,indx+1]], _extra=extra
    endif

return
end    
