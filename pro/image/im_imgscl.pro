;+
; NAME:
;   IM_IMGSCL()
;
; PURPOSE:
;   Scale/stretch an image for display.
;
; INPUTS:
;   image - input image
;
; OPTIONAL INPUTS:
;   losig - low sigma threshold (default -2)
;   hisig - high sigma threshold (default 3)
;   boxfrac - fraction of the image in each dimension to use when
;     computing statistics (default 0.2)
;   topvalue - top image scaling value (default 239) 
;   minvalue - bottom image scaling value (default -10)
;
; KEYWORD PARAMETERS:
;   log - take the log of the image
;   sqrroot - take the square root of the image
;   negative - invert the image
;
; OUTPUTS:
;   img - byte-scaled image
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Mar 22, U of A - written
;   jm05jul25uofa - several improvements
;   jm11aug08ucsd - documented
;
; Copyright (C) 2005, 2011, John Moustakas
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

function im_imgscl, image, losig=losig, hisig=hisig, boxfrac=boxfrac, $
  log=log, sqrroot=sqrroot, negative=negative, topvalue=topvalue, $
  minvalue=minvalue, sigrej=sigrej
    
    imsize = size(image,/dimension)
    xsize = imsize[0] & xcen = xsize/2.0
    ysize = imsize[1] & ycen = ysize/2.0
    
    if (n_elements(losig) eq 0L) then losig = -2.0
    if (n_elements(hisig) eq 0L) then hisig = 3.0
    if (n_elements(boxfrac) eq 0L) then boxfrac = 0.20
    if (n_elements(topvalue) eq 0L) then topvalue = 239L ; !d.table_size-2
    if (n_elements(minvalue) eq 0L) then minvalue = -10L
    if (n_elements(sigrej) eq 0) then sigrej = 5.0
    
    xbox = fix(xsize*boxfrac)/2L
    ybox = fix(ysize*boxfrac)/2L

    im = image
    if keyword_set(log) then im = alog10(float(image))
    if keyword_set(sqrroot) then im = sqrt(float(image>0))

    stats = im_stats(im,sigrej=sigrej)
    substats = im_stats(im[xcen-xbox+1L:xcen+xbox-1L,ycen-ybox+1L:ycen+ybox-1L],sigrej=sigrej)

    mmin = stats.minrej
    mmax = stats.maxrej

;   mmin = (mn+losig*rms);>min(im)
;   mmax = (mn+hisig*rms);<max(im)
    
    img = imgscl(im,min=mmin,max=mmax,top=topvalue)
    if keyword_set(negative) then img = bytscl(topvalue-img,$
      min=minvalue,max=topvalue)

; Hogg's image scaling
    
;   img = ((lo*rms+mean)-float(image))/(lo*rms-hi*rms) ; negative image
;   img = byte((floor(img*255.99) > 0) < 255)

; Christy's image scaling    
    
;   topvalue = 250L
;   minvalue = -40L

;   img = imgscl(image,min=(mean+lo*rms)>min(image),max=(mean+hi*rms)<max(image),top=topvalue)
;   img = bytscl(image,min=min(image),max=max(image),top=topvalue)
;   img = bytscl(topvalue-img,min=minvalue,top=topvalue)

;   img = imgscl(image,min=(mean+lo*rms)>(min(image)*0.5),max=(mean+hi*rms)<max(image))

; ATV's image scaling    
    
;   imgmin = min(image[xcen-xbox:xcen+xbox,ycen-ybox:ycen+ybox])
;   imgmax = max(image[xcen-xbox:xcen+xbox,ycen-ybox:ycen+ybox])
;   imgmin = min(image) & imgmax = max(image)
;   offset = imgmin - (imgmax-imgmin)*0.01

;   img = bytscl(alog10(image-offset),min=alog10(imgmin-offset),$
;     max=alog10(imgmax-offset),top=!d.table_size-2,/nan)

return, img
end    

