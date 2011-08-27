;+
; NAME:
;   IM_HESSPLOT
;
; PURPOSE:
;   Generate a Hess diagram.
;
; INPUTS:
;   x - X-values
;   y - Y-values
;
; OPTIONAL INPUTS:
;   logexp    - logarithmic scaling factor [default 1.0; see
;               LOGSCL()]
;   nxbin2d   - number of X-pixels in the output image
;   nybin2d   - number of Y-pixels in the output image
;   minpts    - if the number of data points in a pixel is less
;               than MINPTS then make the pixel value zero
;               (default 1)
;   extra     - extra keywords for PLOTIMAGE
;
; KEYWORD PARAMETERS:
;   negative  - generate a negative image (e.g., for postscript
;               output) 
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;   image - binned image [NXBIN2D,NYBIN2D]
;
; COMMENTS:
;
; EXAMPLES:
;   IDL> x = randomn(seed,1E5) & y = randomn(seed,1E5)
;   IDL> im_hessplot, x, y, xrange=[-2.5,2.5], yrange=[-4,4], /negative
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 July 24, U of A - written
;   jm08jun25nyu - allow rectangular images; choose the default
;                  dimensions as in HOGG_SCATTERPLOT
;
; Copyright (C) 2006, 2008, John Moustakas
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

pro im_hessplot, x, y, image, logexp=logexp, nxbin2d=nxbin2d, $
  nybin2d=nybin2d, minpts=minpts, negative=negative, _extra=extra

    npts = n_elements(x)
    if (npts eq 0L) then begin
       doc_library, 'im_hessplot'
       return
    endif

    if (n_elements(logexp) eq 0L) then logexp = 1.0
    if (n_elements(nxbin2d) eq 0L) then nxbin2d = ceil(0.3*sqrt(npts))>10L
    if (n_elements(nybin2d) eq 0L) then nybin2d = ceil(0.3*sqrt(npts))>10L
    if (n_elements(minpts) eq 0L) then minpts = 1.0
    
    if size(extra,/type) ne 8L then extra = {dummy: ''}
    
    if tag_exist(extra,'CHARSIZE') then charsize = extra.charsize else charsize = 2.0
    if tag_exist(extra,'XSTYLE') then xstyle = extra.xstyle else xstyle = 3
    if tag_exist(extra,'YSTYLE') then ystyle = extra.ystyle else ystyle = 3
    if tag_exist(extra,'XRANGE') then xrange = extra.xrange else xrange = minmax(x) ; [0.0,1.0]
    if tag_exist(extra,'YRANGE') then yrange = extra.yrange else yrange = minmax(y) ; [0.0,1.0]
    if tag_exist(extra,'XTITLE') then begin
       xtitle = textoidl(extra.xtitle)
       extra = struct_trimtags(extra,except='XTITLE')
    endif else xtitle = ''
    if tag_exist(extra,'YTITLE') then begin
       ytitle = textoidl(extra.ytitle) 
       extra = struct_trimtags(extra,except='YTITLE')
    endif else ytitle = ''

; populate the points in a 2D image
    xbin = abs(float(xrange[1])-xrange[0]) / nxbin2d
    ybin = abs(float(yrange[1])-yrange[0]) / nybin2d

    if (xrange[1] lt xrange[0]) then $
      ximage = (xrange[0] - x) / xbin else $
        ximage = (x - xrange[0]) / xbin
    if (yrange[1] lt yrange[0]) then $
      yimage = (yrange[0] - y) / ybin else $
        yimage = (y - yrange[0]) / ybin

    nx = fix(abs(xrange[1]-xrange[0])/xbin)
    ny = fix(abs(yrange[1]-yrange[0])/ybin)

    image = fltarr(nx,ny)
    populate_image, image, ximage, yimage, assign='cic'

    toofew = where((image gt 0.0) and (image le minpts),ntoofew)
    if (ntoofew ne 0L) then image[toofew] = 0.0

; load the appropriate color table; byte-scale the image
; logarithmically, setting the background to be white (color-table
; code thanks to C. Tremonti)
    loadct, 0, /silent
    tvlct, r, g, b, /get
    cmin = 240 & cmax = 255
    rr = r & gg = g & bb = b
    rr[cmin:cmax] = 255B & gg[cmin:cmax] = 255B & bb[cmin:cmax] = 255B
    tvlct, rr, gg, bb

    img = logscl(image,mean=0.5,exponent=logexp,negative=negative)
    white = where(img eq 255B,nwhite)
    if (nwhite ne 0L) then img[white] = 240B
    
    plotimage, img, charsize=charsize, xstyle=xstyle, ystyle=ystyle, $
      imgxrange=xrange, imgyrange=yrange, xtitle=xtitle, ytitle=ytitle, $
      _extra=extra

return
end
    
