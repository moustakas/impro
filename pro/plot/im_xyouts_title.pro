;+
; NAME:
;   IM_XYOUTS_TITLE
;
; PURPOSE:
;   Flexibly render a plot title on an existing plot using XYOUTS. 
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   xtitle - x-title to use
;   ytitle - y-title to use
;   title - title to use
;   charsize - character size (default !p.charsize>1)
;   xspacing - multiplicative scale factor for YTITLE
;   yspacing - multiplicative scale factor for XTITLE
;   tspacing - multiplicative scale factor for TITLE
;   extra - extra keywords for XYOUTS
;
; KEYWORD PARAMETERS: 
;   xlog - existing plot is logarithmically spaced in x
;   ylog - existing plot is logarithmically spaced in y
;   
; OUTPUTS: 
;   Flexibly adds a title to an existing plot. 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Apr 21, U of A
;   jm07mar15nyu - added [X,Y]LOG keywords
;   jm11sep05ucsd - TITLE optional input added
;
; Copyright (C) 2005, 2007, 2011, John Moustakas
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

pro im_xyouts_title, xtitle=xtitle, ytitle=ytitle, title=title, charsize=charsize, $
  xspacing=xspacing, yspacing=yspacing, xlog=xlog, ylog=ylog, _extra=extra

    if (n_elements(charsize) eq 0) then charsize = !p.charsize>1
    if (n_elements(xspacing) eq 0) then xspacing = 8.0
    if (n_elements(yspacing) eq 0) then yspacing = 3.0
    if (n_elements(tspacing) eq 0) then tspacing = 1.1
    
    if keyword_set(xlog) then xrange = 10^!x.crange else xrange = !x.crange
    if keyword_set(ylog) then yrange = 10^!y.crange else yrange = !y.crange

    if (n_elements(xtitle) ne 0) then begin
       xypos = convert_coord(mean(xrange),yrange[0],/data,/to_normal)
       xpos = reform(xypos[0,*])
       ypos = reform(xypos[1,*]) + (!d.y_ch_size/float(!d.y_size))*(yspacing>charsize)

       xyouts, xpos, ypos, textoidl(xtitle), /normal, align=0.5, $
         charsize=charsize, _extra=extra
    endif

    if (n_elements(ytitle) ne 0) then begin
       xypos = convert_coord(xrange[1],mean(yrange),/data,/to_normal)
       xpos = reform(xypos[0,*]) + (!d.x_ch_size/float(!d.x_size))*(xspacing>charsize)
       ypos = reform(xypos[1,*]) 

       xyouts, xpos, ypos, textoidl(ytitle), /normal, align=0.5, $
         orientation=90, charsize=charsize, _extra=extra
    endif

    if (n_elements(title) ne 0) then begin
       xypos = convert_coord(mean(xrange),yrange[1],/data,/to_normal)
       xpos = reform(xypos[0,*])
       ypos = reform(xypos[1,*]) + (!d.y_ch_size/float(!d.y_size))*(tspacing>charsize)

       xyouts, xpos, ypos, textoidl(title), /normal, align=0.5, $
         charsize=charsize, _extra=extra
    endif

return
end
    
    
