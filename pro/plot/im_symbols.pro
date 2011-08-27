;+
; NAME:
;       IM_SYMBOLS
;
; PURPOSE:
;       Define a variety of user plot symbols.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       psym  - plot symbol number
;
; OPTIONAL INPUTS:
;       psize - symbol size (default 1.0)
;       thick - symbol thickness (default 1.0)
;       color - symbol color as defined in DJS_ICOLOR() (default none)  
;
; KEYWORD PARAMETERS:
;       fill - fill the user symbol
;       help - list and plot all the available symbols 
;
; OUTPUTS:
;       The appropriate plot symbol is stored in PSYM=8 (see the 
;       documentation for PLOTSYM).
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       USERSYM, DJS_ICOLOR()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 14, U of A - written 
;       jm05jun21uofa - documented
;
; Copyright (C) 2002, 2005, John Moustakas
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
;
;-

pro im_symbols, psym, psize=psize, thick=thick, color=color, fill=fill, help=help

    symname = [ $
      'plus_sign', $            ; 101
      'asterisk', $             ; 102
      'point', $                ; 103
      'diamond', $              ; 104
      'triangle', $             ; 105
      'square', $               ; 106
      'X', $                    ; 107
      'circle', $               ; 108
      'triangle_180', $         ; 109
      'pointed_star', $         ; 110
      'down_arrow', $           ; 111
      'up_arrow', $             ; 112
      'left_arrow', $           ; 113
      'right_arrow', $          ; 114
      'star', $                 ; 115
      'circle_0', $             ; 116
      'circle_45', $            ; 117
      'circle_90', $            ; 118
      'circle_180', $           ; 119
      'circle_270', $           ; 120
      'rectangle', $            ; 121
      'cross', $                ; 122
      'xbox', $                 ; 123
      'pentagon', $             ; 124
      'hexagon' $               ; 125
      ]
    nsym = n_elements(symname)
    symindx = lindgen(nsym)
    psym2sym = symindx + 101L   ; symbol number

    if keyword_set(help) then begin

       print, 'Syntax - im_symbols, psym, [psize=,thick=,fill=,color=,help=help]'
       
       len = strtrim(max(strlen(symname)))
;      for i = 0L, nsym-1L do print, psym2sym[i], symname[i], format='(I3,1x,A'+len+')'

       half = ceil(nsym/2.0)
       
       window, xsize=600, ysize=40.*half, /free

       plot, [0], [0], /nodata, xmargin=[0,0], ymargin=[0,0], xsty=5, ysty=5

       for i = 0, half-1L do begin
          im_symbols, psym2sym[i], psize=2.0
          plots, 0.05, 1.0-(i+1.0)/(half+1), ps=8, /norm
          im_symbols, psym2sym[i], psize=2.0, /fill
          plots, 0.1, 1.0-(i+1.0)/(half+1), ps=8, /norm
          xyouts, 0.15, 1.0-(i+1.0)/(half+1), strupcase(symname[i])+' - '+string(psym2sym[i],format='(I0)'), $
            align=0.0, /norm, charthick=2.0, charsize=1.3
       endfor

       i = 0
       for j = half, nsym-1 do begin
          im_symbols, psym2sym[j], psize=2.0
          plots, 0.6, 1.0-(i+1.0)/(half+1), ps=8, /norm
          im_symbols, psym2sym[j], psize=2.0, /fill
          plots, 0.65, 1.0-(i+1.0)/(half+1), ps=8, /norm
          xyouts, 0.7, 1.0-(i+1.0)/(half+1), strupcase(symname[j])+' - '+string(psym2sym[j],format='(I0)'), $
            align=0.0, /norm, charthick=2.0, charsize=1.3
          i = i + 1
       endfor
       
       retall

    endif
    
    npsym = n_elements(psym)
    if npsym eq 0L then im_symbols, /help else psym = psym[0]
    if n_elements(psize) eq 0L then psize = 1.0

    if size(psym,/type) eq 7 then begin

       case psym of
          'plus_sign':  begin   ; plus sign
             yarr = [0, 0, 0, -1, 1] * psize[0]
             xarr = [1, -1, 0, 0, 0] * psize[0]
             fill = 0
          end
          'asterisk':  begin    ; asterisk
             yarr = [-1, 1, 0, 1, -1,0,1,-1,0,0,0] * psize[0]
             xarr = [-1, 1, 0, -1, 1,0,0,0,0,1,-1] * psize[0]
             fill = 0
          end
          'point':  begin       ; point
             yarr = [0,0] * psize[0]
             xarr = [0,0] * psize[0]
             fill = 0
          end
          'diamond':  begin     ; diamond
             yarr = [0, -1, 0, 1, 0] * psize[0]
             xarr = [1, 0, -1, 0, 1] * psize[0]
          end
          'triangle':  begin    ; triangle
             xarr = [-1,0,1,-1]*psize[0]
             yarr = [-1,1,-1,-1]*psize[0]
          end
          'square':  begin      ; square
             yarr = [1, -1, -1,  1, 1] * psize[0]
             xarr = [1,  1, -1, -1, 1] * psize[0]
          end
          'X': begin            ; X
             yarr = [-1, 1, 0, -1, 1] * psize[0]
             xarr = [1, -1, 0, -1, 1] * psize[0]
          end
          'circle':  begin      ; circle
             ang = 2*!PI*FINDGEN(49)/48. ; get position every 5 deg
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          'triangle_180':  begin ; upside down triangle
             xarr = [-1, 0, 1, -1]*psize[0]
             yarr = [ 1,-1, 1, 1]*psize[0]
          end
          'pointed_star':  begin ; pointed star
             r = psize[0]
             ang = (720. / 5*FINDGEN(6) + 45) / !RADEG ; define star angles every 144 deg
             xarr = r*COS(ang)
             yarr = r*SIN(ang)
          end
          'down_arrow': begin   ; down arrow
             xarr = [0,0,.5,0,-.5]*psize[0]
             yarr = [0,-2,-1.4,-2,-1.4]*psize[0]
             fill = 0
          end
          'up_arrow': begin     ; up arrow
             xarr = [0,0,.5,0,-.5]*psize[0]
             yarr = [0,2,1.4,2,1.4]*psize[0]
             fill = 0
          end
          'left_arrow': begin   ; left pointing arrow
             yarr = [0, 0, 0.5, 0, -.5]*psize[0]
             xarr = [0,-2,-1.4,-2,-1.4]*psize[0]
             fill=0
          end
          'right_arrow': begin  ; right pointing arrow
             yarr = [0, 0, 0.5, 0, -.5] * psize[0]
             xarr = [0, 2, 1.4, 2, 1.4] * psize[0]
             fill=0
          end
          'star':  begin        ; star without inlines
             yarr=[1,0.32,0.32,-0.12,-0.81,-0.32,-0.81,-0.12,0.32,0.32,1] *psize[0]
             xarr=[0,0.25,0.94,0.37,0.59,0,-0.59,-0.37,-0.94,-0.25,0] *psize[0]
          end
          'circle_0':  begin    ; circle half filled 0°
             ang=FINDGEN(31)*12.+0.
             ang=[ang,0.+180.]
             ang=[ang,0.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)  &  yarr = psize[0]*SIN(ang)
          end
          'circle_45':  begin   ; circle half filled 45°
             ang=FINDGEN(31)*12.+45.
             ang=[ang,45.+180.]
             ang=[ang,45.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          'circle_90':  begin   ; circle half filled 90°
             ang=FINDGEN(31)*12.+90.
             ang=[ang,90.+180.]
             ang=[ang,90.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          'circle_180':  begin  ; circle half filled 180°
             ang=FINDGEN(31)*12.+180.
             ang=[ang,180.+180.]
             ang=[ang,180.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          'circle_270':  begin  ; circle half filled 270°
             ang=FINDGEN(31)*12.+270.
             ang=[ang,270.+180.]
             ang=[ang,270.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          'rectangle':  begin  ;Square half filled
             xarr = [-1,-1, 1, 1,-1,-1, 1, 1,-1,-1] * psize[0]
             yarr = [ 0, 1, 1, 0, 0,-1,-1, 1, 1, 0] * psize[0]
          end
          'cross': begin
             xarr = [-1,-.3,-.3,.3,.3,1,1,.3,.3,-.3,-.3,-1,-1] * psize[0]
             yarr = [.3,.3,1,1,.3,.3,-.3,-.3,-1,-1,-.3,-.3,.3] * psize[0]
          end
          'xbox': begin
             xarr = [-1,-1,1,-1,1,-1,1,1] * psize[0]
             yarr = [-1,1,1,-1,-1,1,1,-1] * psize[0]
          end
          'pentagon': begin
             xarr = (2*[.5,1,.85,.15,0,.5]-1) * psize[0]
             yarr = (2*[1,.6,0,0,.6,1]-1) * psize[0]
          end
          'hexagon': begin
             xarr = (2*[1,.75,.25,0,.25,.75,1]-1) * psize[0]
             yarr = (2*[.5,1,1,.5,0,0,.5]-1) * psize[0]
          end

          else: begin ; default circle
             ang = 2*!PI*FINDGEN(49)/48. ; get position every 5 deg
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          endelse
          
       endcase 
       
    endif else begin

       psym = long(psym)

; check for and re-assign pre-defined IDL psym keywords

       case psym of
          1L: psym = 101L
          2L: psym = 102L
          3L: psym = 103L
          4L: psym = 104L
          5L: psym = 105L
          6L: psym = 106L
          7L: psym = 107L
          else:
       endcase

       case psym of
          101L:  begin   ; plus sign
             yarr = [0, 0, 0, -1, 1] * psize[0]
             xarr = [1, -1, 0, 0, 0] * psize[0]
             fill = 0
          end
          102L:  begin    ; asterisk
             yarr = [-1, 1, 0, 1, -1,0,1,-1,0,0,0] * psize[0]
             xarr = [-1, 1, 0, -1, 1,0,0,0,0,1,-1] * psize[0]
             fill = 0
          end
          103L:  begin       ; point
             yarr = [0,0] * psize[0]
             xarr = [0,0] * psize[0]
             fill = 0
          end
          104L:  begin     ; diamond
             yarr = [0, -1, 0, 1, 0] * psize[0]
             xarr = [1, 0, -1, 0, 1] * psize[0]
          end
          105L:  begin    ; triangle
             xarr = [-1,0,1,-1]*psize[0]
             yarr = [-1,1,-1,-1]*psize[0]
          end
          106L:  begin        ; square
             yarr = [1, -1, -1,  1, 1] * psize[0]
             xarr = [1,  1, -1, -1, 1] * psize[0]
          end
          107L: begin            ; X
             yarr = [-1, 1, 0, -1, 1] * psize[0]
             xarr = [1, -1, 0, -1, 1] * psize[0]
          end
          108L:  begin      ; circle
             ang = 2*!PI*FINDGEN(49)/48. ; get position every 5 deg
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          109L:  begin ; upside down triangle
             xarr = [-1, 0, 1, -1]*psize[0]
             yarr = [ 1,-1, 1, 1]*psize[0]
          end
          110L:  begin ; pointed star
             r = psize[0]
             ang = (720. / 5*FINDGEN(6) + 45) / !RADEG ; define star angles every 144 deg
             xarr = r*COS(ang)
             yarr = r*SIN(ang)
          end
          111L: begin   ; down arrow
             xarr = [0,0,.5,0,-.5]*psize[0]
             yarr = [0,-2,-1.4,-2,-1.4]*psize[0]
             fill = 0
          end
          112L: begin     ; up arrow
             xarr = [0,0,.5,0,-.5]*psize[0]
             yarr = [0,2,1.4,2,1.4]*psize[0]
             fill = 0
          end
          113L: begin   ; left pointing arrow
             yarr = [0, 0, 0.5, 0, -.5]*psize[0]
             xarr = [0,-2,-1.4,-2,-1.4]*psize[0]
             fill=0
          end
          114L: begin  ; right pointing arrow
             yarr = [0, 0, 0.5, 0, -.5] * psize[0]
             xarr = [0, 2, 1.4, 2, 1.4] * psize[0]
             fill=0
          end
          115L:  begin        ; star without inlines
             yarr=[1,0.32,0.32,-0.12,-0.81,-0.32,-0.81,-0.12,0.32,0.32,1] *psize[0]
             xarr=[0,0.25,0.94,0.37,0.59,0,-0.59,-0.37,-0.94,-0.25,0] *psize[0]
          end
          116L:  begin    ; circle half filled 0°
             ang=FINDGEN(31)*12.+0.
             ang=[ang,0.+180.]
             ang=[ang,0.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)  &  yarr = psize[0]*SIN(ang)
          end
          117L:  begin   ; circle half filled 45°
             ang=FINDGEN(31)*12.+45.
             ang=[ang,45.+180.]
             ang=[ang,45.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          118L:  begin   ; circle half filled 90°
             ang=FINDGEN(31)*12.+90.
             ang=[ang,90.+180.]
             ang=[ang,90.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          119L:  begin  ; circle half filled 180°
             ang=FINDGEN(31)*12.+180.
             ang=[ang,180.+180.]
             ang=[ang,180.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          120L:  begin  ; circle half filled 270°
             ang=FINDGEN(31)*12.+270.
             ang=[ang,270.+180.]
             ang=[ang,270.+180.-FINDGEN(15)*12.]
             ang=ang/180.*!PI
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          end
          121L:  begin  ;Square half filled
             xarr = [-1,-1, 1, 1,-1,-1, 1, 1,-1,-1] * psize[0]
             yarr = [ 0, 1, 1, 0, 0,-1,-1, 1, 1, 0] * psize[0]
          end
          122L: begin
             xarr = [-1,-.3,-.3,.3,.3,1,1,.3,.3,-.3,-.3,-1,-1] * psize[0]
             yarr = [.3,.3,1,1,.3,.3,-.3,-.3,-1,-1,-.3,-.3,.3] * psize[0]
          end
          123L: begin
             xarr = [-1,-1,1,-1,1,-1,1,1] * psize[0]
             yarr = [-1,1,1,-1,-1,1,1,-1] * psize[0]
          end
          124L: begin
             xarr = (2*[.5,1,.85,.15,0,.5]-1) * psize[0]
             yarr = (2*[1,.6,0,0,.6,1]-1) * psize[0]
          end
          125L: begin
             xarr = (2*[1,.75,.25,0,.25,.75,1]-1) * psize[0]
             yarr = (2*[.5,1,1,.5,0,0,.5]-1) * psize[0]
          end

          else: begin ; default circle
             ang = 2*!PI*FINDGEN(49)/48. ; get position every 5 deg
             xarr = psize[0]*COS(ang)
             yarr = psize[0]*SIN(ang)
          endelse
          
       endcase 

    endelse

    colorsym = djs_icolor(color)
    usersym, xarr, yarr, thick=thick, fill=fill, color=colorsym
    
return
end    
