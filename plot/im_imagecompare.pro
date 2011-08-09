;+
; NAME:
;	IM_IMAGECOMPARE
;
; PURPOSE:
;	Interactively compare two images in a single graphics window.
;
; INPUTS:
;	image1 : 2D image
;	image2 : second 2D image to compare to the first
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Use any of the mouse buttons to move the second image around,
;	and press the first and third mouse buttons together to quit.
;
; COMMENTS:
;	Any sized images can be compared, but the larger image needs
;	to be larger in both dimensions.  The images can be of any
;	datatype.  
;
; PROCEDURES USED:
;	DATATYPE()	
;
; MODIFICATION HISTORY:
;	John Moustakas, August 2000, UCB/UofA (with contributions by
;	  Doug Finkbeiner, UCB) 
;
; Copyright (C) 2000, John Moustakas
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

pro im_imagecompare, image1, image2

	on_error, 2	; return to user

	if (not keyword_set(image1)) or (not keyword_set(image2)) then $
          message, 'Two images must be loaded!'

        ndim1 = size(image1,/n_dim)
        ndim2 = size(image2,/n_dim)

        if (ndim1 or ndim2) ne 2 then $
          message, 'Both images must be 2-dimensional!'

; check the sizes

        imsz1 = size(image1)
        imsz2 = size(image2)

        if ((imsz1[1] gt imsz2[1]) and (imsz1[2] lt imsz2[2])) or $
          ((imsz1[2] gt imsz2[2]) and (imsz1[1] lt imsz2[1])) then $
          message, 'One image needs to be larger in both dimensions!'
        if ((imsz2[1] gt imsz1[1]) and (imsz2[2] lt imsz1[2])) or $
          ((imsz2[2] gt imsz1[2]) and (imsz2[1] lt imsz1[1])) then $
          message, 'One image needs to be larger in both dimensions!'
        
        if (imsz1[1] gt imsz2[1]) and (imsz1[2] gt imsz2[2]) then begin
            bimage = image1	; larger image
            simage = image2	; smaller image
            bxysize = [imsz1[1],imsz1[2]]
            sxysize = [imsz2[1],imsz2[2]]
        endif else begin
            bimage = image2	; larger image
            simage = image1	; smaller image
            bxysize = [imsz2[1],imsz2[2]]
            sxysize = [imsz1[1],imsz1[2]]
        endelse 

; some useful quantities

        xydiff = (bxysize-sxysize)/2L
        sxhalf = sxysize[0]/2L
        syhalf = sxysize[1]/2L
        xyrem = sxysize mod 2L
        
; if the images are not type byte, then scale the colortable relative
; to the larger image.  check for complex datatype.
        
        imin = median(bimage)-2.*stdev(bimage) > min(bimage)
        imax = median(bimage)+10.*stdev(bimage) < max(bimage)
        ncolors = !d.table_size

        if datatype(bimage) ne 'BYT' then begin ; byte
            if datatype(bimage) eq 'COM' then $ ; complex
              bim = bytscl(abs(bimage),min=imin,max=imax,top=ncolors-1L) else $
              bim = bytscl(bimage,min=imin,max=imax,top=ncolors-1L)
        endif else bim = bimage
        if datatype(simage) ne 'BYT' then begin
            if datatype(simage) eq 'COM' then $
              sim = bytscl(abs(simage),min=imin,max=imax,top=ncolors-1L) else $
              sim = bytscl(simage,min=imin,max=imax,top=ncolors-1L)
        endif else sim = simage 

; embed the smaller image in a blank image the size of the larger image
        
        eim = bytarr(bxysize[0],bxysize[1])
        eim[xydiff[0]:xydiff[0]+sxysize[0]-1L,$
            xydiff[1]:xydiff[1]+sxysize[1]-1L] = sim 

; display the images

;       window, 2, xsize=(300>bxysize[0])<600, ysize=(300>bxysize[1])<600
        window, 2, xsize=bxysize[0], ysize=bxysize[1]
        tv, bim + eim
        
        print
        print, 'Press the first and third mouse buttons down to quit: '
        print 
        while (1L eq 1L) do begin

            cursor, x, y, 1, /device
            
; check x-coordinates
            
            if (x-sxhalf) lt 0L then x1 = abs(x-sxhalf) else x1 = 0L
            if (x+sxhalf) gt bxysize[0]-1L then $
              x2 = sxysize[0]-(x+sxhalf+xyrem[0]-bxysize[0])-1L else $
              x2 = 2L*sxhalf-1L

; check y-coordinates
            
            if (y-syhalf) lt 0L then y1 = abs(y-syhalf) else y1 = 0L
            if (y+syhalf) gt bxysize[1] then $
              y2 = sxysize[1]-(y+syhalf+xyrem[1]-bxysize[1])-1L else $
              y2 = 2*syhalf-1L

            xstart = (x-sxhalf) > 0L & xend = (x+sxhalf) < bxysize[0] - 1L
            ystart = (y-syhalf) > 0L & yend = (y+syhalf) < bxysize[1] - 1L

            simshift = bytarr(bxysize[0],bxysize[1]) ; shifted smaller image
            simshift[xstart:xend,ystart:yend] = sim[x1:x2,y1:y2]

            tv, bim + simshift ; redisplay

            print, x, y, !err
            if (!err and 5) eq 5 then return

        endwhile

return
end








