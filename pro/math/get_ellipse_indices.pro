;+
; NAME:
;   GET_ELLIPSE_INDICES()
;
; PURPOSE:
;   Identify the elements lying inside a given ellipse.
;
; INPUTS: 
;   x, y - input values
;
; OPTIONAL INPUTS: 
;   major, minor - axis lengths in the same units as X,Y (default 1.0) 
;   angle - angle of the ellipse, measured counter-clockwise from the
;     X axis (default 0)
;   xcenter, ycenter - origin of the ellipse (default 0)
;
; KEYWORD PARAMETERS: 
;   debug - make a simple debugging plot
;
; OUTPUTS: 
;   indx - indices of X,Y values lying inside the specified ellipse 
;
; COMMENTS:
;   This routine can be used in tandem with the output of
;   COVAR2ELLIPSE(). 
;
;   The angle is defined as in TVELLIPSE.  Portions of the algorithm
;   were taken from DIST_ELLIPSE.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jul 09, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function get_ellipse_indices, x, y, major=major, minor=minor, $
  angle=angle, xcenter=xcenter, ycenter=ycenter, debug=debug

    npts = n_elements(x)
    if (npts eq 0) or (npts ne n_elements(y)) then begin
       doc_library, 'get_ellipse_indices'
       return, -1
    endif
    
    if (n_elements(xcenter) eq 0) then xcenter = 0.0
    if (n_elements(ycenter) eq 0) then ycenter = 0.0
    
; rotate pixels to match ellipse orientation
    ang = (90-angle)*!dtor
    rot_ccw = [ [cos(ang),sin(ang)], [-sin(ang),cos(ang)] ]
    xyprime = transpose([[(x-xcenter)],[(y-ycenter)]]) ## rot_ccw

    indx = where(sqrt((xyprime[0,*]*major/minor)^2 + xyprime[1,*]^2) lt major,nindx)

    if keyword_set(debug) then begin
       djs_plot, xyprime[0,*], xyprime[1,*], ps=6
       if (nindx ne 0) then djs_oplot, xyprime[0,indx], $
         xyprime[1,indx], ps=6, color='orange'
    endif

return, indx
end
