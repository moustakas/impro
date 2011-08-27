;+
; NAME:
;       IM_OPENCLOSE
;
; PURPOSE:
;       Simple wrapper for opening and closing a postscript file
;       (calls DFPSPLOT/DFPSCLOSE), or for waiting for a keystroke. 
;
; INPUTS:
;       psname - postscript file name (without .ps or .eps extension)
;
; OPTIONAL INPUTS:
;       extra - keywords for DFPSPLOT
;
; KEYWORD PARAMETERS:
;       square       - square postscript output
;       postscript   - open or close the postscript file 
;       encapsulated - generate an encapsulated postscript file and
;                      append a .eps suffix, otherwise the output
;                      suffix is set to .ps
;       png          - NOT SUPPORTED YET
;       close        - close the postscript file
;       silent       - do not print to STDOUT
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The PNG keyword is stronger than the postscript keyword. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, longtimeago, U of A - written
;       jm04mar10uofa - do not ISOLATIN1 fonts; remove !p.font calls 
;       jm05mar30uofa - added ENCAPSULATED keyword; figure out the
;                       postscript file extension automatically
;       jm05may18uofa - added PNG keyword
;
; Copyright (C) 2004-2005, John Moustakas
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

pro im_openclose, psname, square=square, postscript=postscript1, $
  encapsulated=encapsulated, png=png, close=close, silent=silent, $
  color=color, _extra=extra
    
    if (n_elements(psname) eq 0L) and (not keyword_set(close)) then begin
       doc_library, 'im_openclose'
       return
    endif

    if (n_elements(color) eq 0L) then color = 1L
    
    if keyword_set(png) then postscript1 = 0L

    if n_elements(square) eq 0L then square = 1L
    if n_elements(extra) ne 0L then if tag_exist(extra,'landscape') then square = 0
    
    if keyword_set(postscript1) then begin

       if keyword_set(encapsulated) then suffix = '.eps' else suffix = '.ps'

       if keyword_set(close) then dfpsclose else begin
          if (strmatch(strtrim(psname,2),'*.ps',/fold) eq 1B) or $
            (strmatch(strtrim(psname,2),'*.eps',/fold) eq 1B) then $
            psname = repstr(repstr(psname,'.ps',''),'.eps','')
          dfpsplot, psname+suffix, color=color, square=square, _extra=extra
          if not keyword_set(silent) then splog, 'Writing '+psname+suffix
       endelse

    endif else if keyword_set(close) and (not keyword_set(png)) then begin

       if not keyword_set(silent) then splog, 'Press any key to continue.'
       cc = strupcase(get_kbrd(1))

    endif

    if keyword_set(png) then begin

       suffix = '.png'

       if keyword_set(close) then begin
          img = tvrd()
          tvlct, r, g, b, /get
          write_png, psname+suffix, img, r, g, b
          if (!d.name ne 'X') then set_plot, 'X'
       endif else begin
          if (strmatch(strtrim(psname,2),'*.png',/fold) eq 1B) then $
            psname = repstr(psname,'.png','')
          set_plot, 'Z'       
          if not keyword_set(silent) then splog, 'Writing '+psname+suffix
       endelse

    endif
    
return
end    
