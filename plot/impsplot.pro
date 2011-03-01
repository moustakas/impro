;+
; NAME:
;       IMPSPLOT
;
; PURPOSE:
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2008, John Moustakas
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

pro impsplot, psfile, xsize=xsize, ysize=ysize, pdf=pdf, close=close, $
  position=position, _extra=extra

    if keyword_set(close) then begin
       if (!d.name eq 'X') then begin                                                              
          print, 'Postscript device not open!'
       endif else begin 
          device, /close
          set_plot, 'X'
       endelse
       return
    endif

stop    
    
    if (n_elements(nx) eq 0L) then nx = 1L
    if (n_elements(ny) eq 0L) then nx = 1L
    if (n_elements(xmargin) eq 0L) then xmargin = [1.1,0.4]
    if (n_elements(ymargin) eq 0L) then ymargin = [0.4,1.1]
    if (n_elements(xspace) eq 0L) then xspace = 0.0
    if (n_elements(yspace) eq 0L) then yspace = 0.0
    if (n_elements(xsize) eq 0L) then xpage = 8.5 else xpage = xsize
    if (n_elements(ysize) eq 0L) then ypage = 8.5 else ypage = ysize
    if (n_elements(width) eq 0L) then width = 6.5
    if (n_elements(height) eq 0L) then height = 6.5

    im_plotfaves, /post
    arm_plotconfig, psfile, nx=nx, ny=ny, landscape=landscape, xspace=xspace, $
      yspace=xspace, width=width, height=height, xmargin=xmargin, $
      ymargin=ymargin, xpage=xpage, ypage=ypage, position=position, $
      /normal
;   pagemaker, nx=nx, ny=ny, landscape=landscape, xspace=xspace, $
;     yspace=xspace, width=width, height=height, xmargin=xmargin, $
;     ymargin=ymargin, xpage=xpage, ypage=ypage, position=position, $
;     /normal
;   set_plot, 'PS'
;   device, file=psfile, xsize=xpage, ysize=ypage, $
;     xoff=0, yoff=0.5+(10-ys), /inch, encap=encap, /color, $
;     _extra=extra

    im_plotfaves

return
end
;------------------------------------------------------------------------------
