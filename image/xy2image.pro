;+
; NAME:
;	XY2IMAGE
;
; PURPOSE:
;	Populate an image with the coordinates given by (x,[y]).
;
; CALLING SEQUENCE:
;	xy2image, x, y=y, image, binsize=binsize, weights=weights
;
; INPUTS:
;	x	: x coordinates(s)
;
; OPTIONALS INPUTS:
;	y	: y coordinate(s)
;	weights	: weight(s) to add at each (x,[y]) position
;	binsize	: resolution of each image pixel
;
; KEYWORD PARAMETERS:
;	float   : will return a float scaled image (default is byte)
;    
; OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; COMMENTS:
;	The binsize should reflect the decimal accuracy of x and y.
;	For example, if the x,y coordinates are known to 0.01 units,
;	then the binsize should be set to 0.01.  
;
; MODIFICATION HISTORY:
;	John Moustakas, 22 August 2000, UofA
;	jm01may30uofa - allowed image to be passed with defined dimensions
;
; Copyright (C) 2000-2001, John Moustakas
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

pro xy2image, x, y=y, image, binsize=binsize, weights=weights, float=float

    if n_params() eq 0L then begin
       print, 'Syntax : xy2image, x, y=y, image, binsize=binsize, weights=weights, float=float'
       return
    endif

    npts = n_elements(x)
    if (keyword_set(y)) then if (npts ne n_elements(y)) then $
      message, 'Dimensions of X and Y do not agree'
    if (keyword_set(weights)) then $
      if (npts ne n_elements(weights)) then $
      message, 'Dimensions of X and WEIGHTS do not agree'
    
    if not keyword_set(y) then y = findgen(npts)
    if not keyword_set(binsize) then binsize = 1L
    if keyword_set(float) then begin
       if not keyword_set(weights) then weights = fltarr(npts) + 1. 
    endif else begin
       if not keyword_set(weights) then weights = bytarr(npts) + 1B
    endelse
    
    if n_elements(image) ne 0L then begin ; fixed image size

       xdim = (long(x/binsize)-1L) > 0L
       ydim = (long(y/binsize)-1L) > 0L

    endif else begin  ; variable image size

       xdim = long((x-min(x))/binsize)
       ydim = long((y-min(y))/binsize)
       
       nx = long(max(xdim)+1L)
       ny = long(max(ydim)+1L)
    
       if keyword_set(float) then image = fltarr(nx,ny) else image = bytarr(nx,ny)
    
    endelse
    
    for k = 0L, npts-1L do image[xdim[k],ydim[k]] = image[xdim[k],ydim[k]] + weights[k]
    
    if not keyword_set(float) then image = image gt 0B
    
return
end







