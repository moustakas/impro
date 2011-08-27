;+
; NAME:
;   DISPLAY_DSS
;
; PURPOSE:
;   Bone-simple script to render an image.
;
; INPUTS: 
;   image - input FITS file name *or* input 2D image 
;
; OPTIONAL INPUTS: 
;   heaader - corresponding image header
;   outfile - output file name (default 'dss_object.png')
;   title - image title
;
; KEYWORD PARAMETERS: 
;   write - screen capture the result and write out a PNG file 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only works on one image at a time.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Dec 11, U of A
;   jm11aug08ucsd - documented
;
; Copyright (C) 2002, 2011, John Moustakas
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

pro display_dss, image, header=header, outfile=outfile, title=title, write=write
    
    nimage = n_elements(image)
    if nimage eq 0L then begin
       print, 'Syntax - display_dss, image'
       return
    endif
    
    imtype = size(image,/type)
    if (n_elements(outfile) eq 0) then outfile = 'dss_object.png'
    if (keyword_set(write) eq 0) then window, 0, xs=650, ys=650

    if (imtype eq 7) then begin

       if file_test(image,/regular) eq 1L then begin

          splog, 'Reading '+image+'.'
          im = readfits(image,header,/silent)

          imsize = size(im,/dimension)
          xsize = imsize[0]
          ysize = imsize[1]

          xmid = float(xsize-1.0)/2.0 & ymid = float(ysize-1.0)/2.0
          pixscale = abs(sxpar(header,'CDELT1'))*3600.0 ; pixel scale [arcsec/pixel]
          
          xaxis = (findgen(xsize)-xmid)*pixscale ;/60.0
          yaxis = (findgen(ysize)-ymid)*pixscale ;/60.0

          xtitle = 'x (arcsec)' & ytitle = 'y (arcsec)'
          
       endif else begin

          splog, 'Image '+image+' not found...returning.'
          return

       endelse

    endif else begin

       im = image

       imsize = size(im,/dimension)
       xsize = imsize[0]
       ysize = imsize[1]

       xaxis = findgen(xsize) & yaxis = findgen(ysize)
       xtitle = 'x (pixel)' & ytitle = 'y (pixel)'
       
    endelse 

    djs_iterstat, im, sigma=imsig, median=immed, sigrej=3.0
;   immin = immed-3*imsig
;   imtop = immed+8*imsig
    
    imdisp, imgscl(im,min=immin,top=imtop,/log), /erase, $
      /axis, /negative, yminor=3, position=[0.15,0.05,0.95,1.0], $
      yrange=minmax(yaxis), xrange=minmax(xaxis), charsize=2.0, charthick=2.0, $
      xtitle=xtitle, ytitle=ytitle, xthick=2.0, ythick=2.0, title=title

; overplot the cardinal directions

    if n_elements(header) ne 0L then $
      arrows, header, 0.8*max(xaxis), 0.7*max(yaxis), /data, thick=4.0
    
; grab the image
    
    if keyword_set(write) then begin
       outim = tvrd()
       splog, 'Writing '+outfile
       write_png, outpath+outfile, outim
    endif else begin
;      print, 'Press any key to continue.'
;      cc = strupcase(get_kbrd(1))
    endelse
       
return
end
