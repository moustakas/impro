;+
; NAME:
;   DS9POLYGON_INDICES()
;
; PURPOSE:
;   Convert a ds9 "polygon" regions file into a list of pixel indices
;   enclosed by that polygon; basically a wrapper on POLYFILLV.
;   Handles both image and celestial coordinates.
;
; INPUTS:
;   regionsfile - ds9 style polygon region file (.reg)
;
; OPTIONAL INPUTS:
;   [x,y]size - size of the image from which REGIONSFILE was produced
;     along the x- and y-dimensions [pixels] 
;                      OR
;   header - FITS header corresponding to the image from which
;     REGIONSFILE was produced 
;	
; KEYWORD PARAMETERS:
;   inverse - return the indices of the pixels *outside* the mask 
;
; OUTPUTS:
;   indices - indices of the pixels enclosed by the polygon
;
; OUTPUTS:
;   xvert, yvert - X and Y vertices (pixels)
;
; COMMENTS:
;   If the regions file is in celestial (ra,dec) coordinates then the
;   (astrometrically complete) FITS header is REQUIRED.  
; 
;   If both [X,Y]SIZE and HEADER are given then HEADER takes
;   precedence.  
;
;  The regions file should contain just one polygon.  (Multiple
;  polygons are not supported because each would in principle contain
;  a different number of indices; better to just call this routine
;  repeatedly.) 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   John Moustakas, 2013 May 29, Siena
;
; Copyright (C) 2013, John Moustakas
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

function ds9polygon_indices, regionsfile, xsize=xsize, ysize=ysize, $
  header=header, xvert=xx, yvert=yy, inverse=inverse

    if n_elements(regionsfile) eq 0 then begin
       doc_library, 'ds9polygon_indices'
       return, -1
    endif

    if file_test(regionsfile) eq 0 then begin
       splog, 'DS9 regions file '+regionsfile+' not found!'
       return, -1
    endif

    nxy = n_elements(xsize)+n_elements(ysize)
    nhead = n_elements(header)
    if nxy eq 0 and nhead eq 0 then begin
       splog, 'Must provide either [X,Y]SIZE *or* HEADER (see also COMMENTS)!'
       return, -1
    endif

    isfk5 = 0 ; celestial or image coordinates?
    haspolygon = 0
    
    openr, lun, regionsfile, /get_lun
    str = ''
    while not eof(lun) do begin
       readf, lun, str
       if strmatch(str,'*fk5*') then isfk5 = 1
       if strmid(str,0,8) eq 'polygon(' then begin
          haspolygon = 1
          xy = double(strsplit(strmid(str,8,strlen(str)-9),',',/extract))
          nvert = n_elements(xy)/2
          ind = lindgen(nvert)
          xx1 = xy[ind*2]
          yy1 = xy[ind*2+1]
          if nhead eq 0 then begin ; no header
             xx = xx1-1 ; ds9 is 1-indexed
             yy = yy1-1
          endif else begin ; has header
             extast, header, astr
             xsize = astr.naxis[0]
             ysize = astr.naxis[1]
; if the polygon vertices are in celestial coordinates, get the
; astrometry from the header
             if isfk5 then begin
                ad2xy, xx1, yy1, astr, xx, yy
             endif else begin
; pixel coordinates
                xx = xx1
                yy = yy1
             endelse
          endelse
          indices = polyfillv(xx,yy,xsize,ysize)
          if keyword_set(inverse) then begin
             all = lindgen(float(xsize)*float(ysize))
             remove, indices, all
             indices = temporary(all)
          endif
       endif 
    endwhile
    free_lun, lun

    if haspolygon eq 0 then begin
       splog, 'No polygon region found in file '+regionsfile+'!'
       return, -1
    endif
    
return, indices
end
    
