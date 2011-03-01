;+
; NAME:
;   IM_GETPOSITION()
;
; PURPOSE:
;   Compute the position coordinates of one or more plots (without
;   using !p.multi) under various constraints.  
;
; INPUTS: 
;   None required.
;
; OPTIONAL INPUTS: 
;   nx - number of plots from left to right (default=1)
;   ny - number of plots from top to bottom (default=1)                  
;   xmargin - left/right margin, array of 1 or 2 values (inches,
;     default=[1.25,0.75]) 
;   ymargin - top/bottom margin, array of 1 or 2 values (inches,
;     default=[0.75,0.75]) 
;   width - width of each plot, array of 1 or nx values (inches) 
;   height - height of each plot, array of 1 or ny values (inches) 
;   xspace - horizontal spacing between plots, array of 1 or nx-1
;     values, (inches) 
;   yspace - vertical spacing between plots, array of 1 or ny-1 values
;     (inches) 
;   xpage - length of plot area from left to right (inches,
;     default=8.5) 
;   ypage - length of plot area from top to bottom (inches,
;     default=11.0)  
;
; KEYWORD PARAMETERS: 
;   landscape - rotates orientation of each plot 90 degrees
;
; OUTPUTS: 
;   position - an N by (nx x ny) array of plot coordinates
;     (normalized), the N elements are [x0,y0,x2,y2] where the plot
;     corners are ordered 0-3 counter-clockwise from the lower-left 
;
; COMMENTS:
;   This routine is designed to be used in tandem with IM_PLOTCONFIG. 
;
; EXAMPLES:
;   plot, findgen(10), position=im_getposition()
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009-Nov-03, UCSD - shamelessly hacked
;     version of Andy Marble's ARM_PLOTCONFIG. 
;
; Copyright (C) 2009, John Moustakas
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

function im_getposition, nx=nx, ny=ny, xpage=xpage, ypage=ypage, $
  xmargin=xmargin, ymargin=ymargin, width=width, height=height, $
  xspace=xspace, yspace=yspace, landscape=landscape

; specify default values
    dflt_xpage =  8.5D           ; xpage default (inches)
    dflt_ypage = 11.0D           ; ypage default (inches)
    dflt_xmargin = 1.0D          ; xmargin default (inches)
    dflt_ymargin = 1.0D          ; ymargin default (inches)
    dflt_nx      = 1D            ; nx default
    dflt_ny      = 1D            ; ny default
    dflt_xspace  = 0.75D         ; xspace default (inches)
    dflt_yspace  = 0.75D         ; yspace default (inches)

; error checking and default values
    if (n_elements(xpage) eq 0) then xpage = dflt_xpage
    if (n_elements(ypage) eq 0) then ypage = dflt_ypage
    
    if (n_elements(xpage) gt 1) or (n_elements(ypage) gt 1) then $
      message, 'XPAGE and YPAGE must be scalars'
    
    xp = xpage
    yp = ypage 
    case n_elements(xmargin) of
       0: xm = [1,1]*dflt_xmargin
       1: xm = [1,1]*xmargin
       2: xm = xmargin
       else: message, 'XMARGIN must have at most 2 elements'
    endcase

    case n_elements(ymargin) of 
       0: ym = [1,1]*dflt_ymargin
       1: ym = [1,1]*ymargin
       2: ym = ymargin
       else: message, 'YMARGIN must have at most 2 elements'
    endcase

    if (n_elements(nx) eq 0) then nx = dflt_nx
    if (n_elements(ny) eq 0) then ny = dflt_ny
    if (n_elements(nx) gt 1) or (nx le 0) then $
      message, 'NX must be a positive scalar'
    if (n_elements(ny) gt 1) or (ny le 0) then $
      message, 'NY must be a positive scalar'
    
    if (nx gt 1) then begin 
       if (n_elements(xspace) eq 0) then xsp = $
         [0.0,replicate(dflt_xspace,nx-1)] else $
           if (n_elements(xspace) eq 1) then xsp = $
         [0.0,replicate(xspace,nx-1)] else $
           if (n_elements(xspace) eq nx-1) then xsp = [0.0,xspace] else $
             message, 'XSPACE must have 1 or NX-1 elements'
    endif else xsp = 0.0

    if (ny gt 1) then begin 
       if (n_elements(yspace) eq 0) then ysp = $
         [0.0,replicate(dflt_yspace,ny-1)] else $
           if (n_elements(yspace) eq 1) then ysp = $
         [0.0,replicate(yspace,ny-1)] else $
           if (n_elements(yspace) eq ny-1) then ysp = [0.0,yspace] else $
             message, 'YSPACE must have 1 or NY-1 elements'
    endif else ysp = 0.0

; calculate widths and heights of plot boxes, if not specified
    if (n_elements(width) eq 0) then ww = $
      replicate((xp-total(xm)-total(xsp))/float(nx),nx) else $
        if (n_elements(width) eq 1) then ww = replicate(width,nx) else $
          if (n_elements(width) ne nx) then $
            message, 'WIDTH must have 1 or NX elements' else ww = width
    if (abs(total(ww)+total(xsp)+total(xm)-xp) gt 1D-3) then $
      message, 'WIDTH and XSPACE values exceed usable area'

    if (n_elements(height) eq 0) then hh = $
      replicate((yp-total(ym)-total(ysp))/float(ny),ny) else $
        if (n_elements(height) eq 1) then hh = replicate(height,ny) else $
          if (n_elements(height) ne ny) then $
            message, 'HEIGHT must have 1 or NY elements' else hh = height
    if (abs(total(hh)+total(ysp)+total(ym)-yp) gt 1D-3) then $
      message, 'HEIGHT and YSPACE values exceed usable area'

; calculate corner coordinates for each plot box 
    position = fltarr(4,nx*ny)
    for ii=0, ny-1 do begin 
       for jj=0, nx-1 do begin
          x0 = xm[0] + total(ww[0:jj]) + total(xsp[0:jj]) - ww[jj] 
          x2 = xm[0] + total(ww[0:jj]) + total(xsp[0:jj])
          y0 = yp - ym[0] - total(hh[0:ii]) - total(ysp[0:ii])
          y2 = yp - ym[0] - total(hh[0:ii]) - total(ysp[0:ii]) + hh[ii] 
          position[*,jj+ii*nx] = [x0,y0,x2,y2]
       endfor
    endfor   
          
; convert physical coordinates to normal coordinates
    position[[0,2],*] = position[[0,2],*] / xp
    position[[1,3],*] = position[[1,3],*] / yp

return, position
end

