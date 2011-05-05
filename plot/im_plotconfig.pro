;+
; NAME:
;   IM_PLOTCONFIG
;
; PURPOSE:
;   Customized wrapper on generating postscript output. 
;
; INPUTS: 
;   plotnum - preset plotting number, specifying the number and
;     orientation of the panels
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

pro im_plotconfig, plotnum, position, psfile=psfile, psclose=psclose, $
  xmargin=xmargin, ymargin=ymargin, xspace=xspace, yspace=yspace, $
  width=width, height=height, xpage=xpage, ypage=ypage, keynote=keynote, $
  blackwhite=blackwhite, pcolor=pcolor, gzip=gzip, pdf=pdf, pskeep=pskeep, $
  silent=silent, _extra=extra
; jm08aug21nyu - simple wrapper on ARM_PLOTCONFIG with a set of
;   predefined plots that I frequently use
; jm09mar20nyu - postscript output gets generated if PSFILE is passed,
;   and depending on its extension (.PS vs .EPS), either regular or 
;   encapsulated postscript is built
; jm09may25nyu - added GZIP and PDF keywords; if /PDF *and* /GZIP are
;   set then /PDF wins; note that both keywords require PSFILE

;   common im_plotconfig, defcolor
;   if (n_elements(defcolor) eq 0) then defcolor = !p.color
    
    if (n_elements(plotnum) eq 0) then plotnum = 0
    case plotnum of
       0: begin ; 1x1 portrait
          nx1 = 1 & ny1 = 1 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 7.0 else width1 = width
          if (n_elements(height) eq 0) then height1 = 7.0 else height1 = height
       end
       1: begin ; 2x1 landscape
          nx1 = 2 & ny1 = 1 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 4.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 4.5 else height1 = height
       end
       2: begin ; 2x2 landscape
          nx1 = 2 & ny1 = 2 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 4.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.5*[1,1] else height1 = height
       end
       3: begin ; 3x1 landscape
          nx1 = 3 & ny1 = 1 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.35,1.05] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.2*[1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.4 else height1 = height
       end
       4: begin ; 1x3 portrait
          nx1 = 1 & ny1 = 3 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 7.0 else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.2*[1,1,1] else height1 = height
       end
       5: begin ; 2x2 portrait
          nx1 = 2 & ny1 = 2 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.5*[1,1] else height1 = height
       end
       6: begin ; 1x2 portrait
          nx1 = 1 & ny1 = 2 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 7.0 else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.5*[1,1] else height1 = height
       end
       7: begin ; 2x3 portrait
          nx1 = 2 & ny1 = 3 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.2*[1,1,1] else height1 = height
       end
       8: begin ; 1x1 landscape
          nx1 = 1 & ny1 = 1 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 9.0 else width1 = width
          if (n_elements(height) eq 0) then height1 = 6.0 else height1 = height
       end
       9: begin ; 1x2 landscape
          nx1 = 1 & ny1 = 2 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 7.0 else width1 = width
          if (n_elements(height) eq 0) then height1 = 4.0*[1,1] else height1 = height
       end
       10: begin ; 2x3 landscape
          nx1 = 2 & ny1 = 3 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 4.75*(intarr(nx1)+1) else width1 = width
          if (n_elements(height) eq 0) then height1 = 2.3*(intarr(ny1)+1) else height1 = height
       end
       11: begin ; 1x3 landscape
          nx1 = 1 & ny1 = 3 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.0] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 9.5*(intarr(nx1)+1) else width1 = width
          if (n_elements(height) eq 0) then height1 = 2.4*(intarr(ny1)+1) else height1 = height
       end
       12: begin ; 2x1 portrait
          nx1 = 2 & ny1 = 1 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 7.0 else height1 = height
       end
       13: begin ; 3x1 portrait
          nx1 = 3 & ny1 = 1 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1D,0.4D] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4D,1.1D] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0D else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0D else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5D*[1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 7.0D else height1 = height
       end
       14: begin ; 3x2 portrait
          nx1 = 3 & ny1 = 2 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.2*[1,1] else height1 = height
       end
       15: begin ; 3x3 portrait
          nx1 = 3 & ny1 = 3 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.2*[1,1,1] else height1 = height
       end
       16: begin ; 3x4 portrait
          nx1 = 3 & ny1 = 4 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.0,0.3] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 2.4*[1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 2.4*[1,1,1,1] else height1 = height
       end
       17: begin ; 3x2 landscape
          nx1 = 3 & ny1 = 2 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = replicate(3.0,nx1) else width1 = width
          if (n_elements(height) eq 0) then height1 = replicate(3.5,ny1) else height1 = height
       end
       18: begin ; 4x2 landscape
          nx1 = 4 & ny1 = 2 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.4,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = replicate(2.3,nx1) else width1 = width
          if (n_elements(height) eq 0) then height1 = replicate(3.5,ny1) else height1 = height
       end
       19: begin ; 2x4 portrait
          nx1 = 2 & ny1 = 4 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 2.4*[1,1,1,1] else height1 = height
       end
       20: begin ; 6x3 landscape
          nx1 = 6 & ny1 = 3 & landscape1 = 1
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.35,1.05] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = replicate(1.6,nx1) else width1 = width
          if (n_elements(height) eq 0) then height1 = replicate(2.0,ny1) else height1 = height
       end
       21: begin ; 4x6 portrait
          nx1 = 4 & ny1 = 6 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [0.8,0.3] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,0.8] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 1.85*(intarr(nx1)+1) else width1 = width
          if (n_elements(height) eq 0) then height1 = 1.65*(intarr(ny1)+1) else height1 = height
       end
       22: begin ; 4x3 portrait
          nx1 = 4 & ny1 = 3 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 1.75*[1,1,1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 3.2*[1,1,1] else height1 = height
       end
       23: begin ; 2x6 portrait
          nx1 = 2 & ny1 = 5 & landscape1 = 0
          if (n_elements(xmargin) eq 0) then xmargin1 = [1.1,0.4] else xmargin1 = xmargin
          if (n_elements(ymargin) eq 0) then ymargin1 = [0.3,1.1] else ymargin1 = ymargin
          if (n_elements(xspace) eq 0) then xspace1 = 0.0 else xspace1 = xspace
          if (n_elements(yspace) eq 0) then yspace1 = 0.0 else yspace1 = yspace
          if (n_elements(width) eq 0) then width1 = 3.5*[1,1] else width1 = width
          if (n_elements(height) eq 0) then height1 = 1.9*[1,1,1,1,1] else height1 = height
       end
       else: begin
          splog, 'Plot number not recognized - write one!'
          stop
       end
    endcase

    xpage1 = total(width1,/double)+total(xmargin1,/double)+total(xspace1,/double)
    ypage1 = total(height1,/double)+total(ymargin1,/double)+total(yspace1,/double)

; get the position vector
    position = im_getposition(nx=nx1,ny=ny1,xmargin=xmargin1,$
      ymargin=ymargin1,width=width1,height=height1,xspace=xspace1,$
      yspace=yspace1,xpage=xpage1,ypage=ypage1,landscape=landscape1)

; close the PS file    
    if keyword_set(psclose) then begin
       if (!d.name ne 'X') then begin
          device, /close
          set_plot, 'X'
          im_plotfaves
;         !p.color = defcolor
          if keyword_set(pdf) or keyword_set(keynote) then begin
             if (n_elements(psfile) eq 0) then begin
                splog, 'You must pass PSFILE to generate a PDF!' 
             endif else begin
                pdffile = repstr(repstr(psfile,'.ps','.pdf'),'.ps','.pdf')
                spawn, 'ps2pdf13 '+psfile+' '+pdffile, /sh
                if (keyword_set(pskeep) eq 0) then rmfile, psfile
             endelse
             return
          endif 
          if keyword_set(gzip) and (n_elements(psfile) ne 0) then begin
             spawn, 'gzip -f '+psfile, /sh
             return
          endif
       endif else splog, 'Postscript device not open!'
    endif

; write out
    if (n_elements(psfile) ne 0) and (keyword_set(psclose) eq 0) then begin
       if (keyword_set(silent) eq 0) then begin
          splog, 'Writing '+psfile
       endif
       im_plotfaves, /postscript, keynote=keynote, _extra=extra

       if strmatch(psfile,'*eps*',/fold) then encap = 1 else encap = 0

       set_plot, 'ps'
       if landscape1 then begin
          xoffset = ((8.5-ypage1)/2.0)>0.0
          if encap then yoffset = xpage1 else $
            yoffset = (11.0-(11.0-xpage1)/2.0)>0.0
       endif else begin    
          xoffset = ((8.5-xpage1)/2.0)>0.0
          yoffset = ((11.0-ypage1)/2.0)>0.0
       endelse
       device, file=psfile, color=(keyword_set(blackwhite) eq 0), /inches, $
         xsize=xpage1, ysize=ypage1, xoffset=xoffset, yoffset=yoffset, $
         landscape=landscape1, encapsulated=encap, _extra=extra
    endif

    if (!d.name eq 'X') then !p.color = djs_icolor('white')
    
;   if keyword_set(keynote) then $
;     !p.color = djs_icolor('white') else $
;       !p.color = djs_icolor('default')
;      if (n_elements(pcolor) eq 0) then $
;        pcolor = fsc_color('white',0) else $
;          pcolor = fsc_color(pcolor)
    
return
end
