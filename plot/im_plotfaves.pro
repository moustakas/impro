;+
; NAME:
;   IM_PLOTFAVES
;
; PURPOSE:
;   Load "favorite" settings for generating plots.
;
; OPTIONAL INPUTS:
;   xticklen : tick length of the x-axis (default 0.03)
;   yticklen : tick length of the y-axis (default 0.03)
;   pthick   : thickness of the lines connecting points (default 2.0)
;   charsize : character size (default 1.8)
;   charthick: character thickness (default 2.0)
;   xthick   : thickness of the x-axis (default 2.0)
;   ythick   : thickness of the y-axis (default 2.0)
;	
; KEYWORD PARAMETERS:
;   postscript: 
;   pdf       : 
;   restore   : restore the plot fields to their IDL default values
;   extras    : load some additional handy plotting variables
;   help      : print the code syntax to the screen and return
;
; COMMENTS:
;   This routine was written to be placed at the beginning of a
;   code which intends to make nice plots which will perhaps be
;   written to postscript.  Note:  this routine should be called
;   again at the end of the code with the /restore keyword set.
;
; MODIFICATION HISTORY:
;   John Moustakas, 2000 September 15, U of A
;   jm08aug21nyu - some tweaks; use CLEANPLOT for restoring; added
;     POSTSCRIPT keyword
;   jm08nov03nyu - added PDF optional input
;   jm09feb10nyu - bug fix when passing optional inputs
;   jm09mar29nyu - added KEYNOTE keyword
;
; Copyright (C) 2000, 2008-2009, John Moustakas
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

pro im_plotfaves, xthick=xthick, ythick=ythick, thick=thick, $
  charthick=charthick, charsize=charsize, color=color, $
  restore=restore, postscript=postscript, keynote=keynote, $
  extras=extras, help=help

    if keyword_set(help) then begin
       doc_library, 'im_plotfaves'
       return
    endif

    if (n_elements(xthick) eq 0)    then !x.thick = 2.0 else !x.thick = xthick
    if (n_elements(ythick) eq 0)    then !y.thick = 2.0 else !y.thick = ythick
    if (n_elements(thick) eq 0)     then !p.thick = 2.0 else !p.thick = thick
    if (n_elements(charthick) eq 0) then !p.charthick = 1.5 else !p.charthick = charthick
    if (n_elements(charsize) eq 0)  then !p.charsize = 1.8 else !p.charsize = charsize
;   if (n_elements(color) eq 0)     then !p.color = fsc_color('white') else !p.color = color

    if keyword_set(postscript) then begin
       if (n_elements(xthick) eq 0)    then !x.thick = 4.0 else !x.thick = xthick
       if (n_elements(ythick) eq 0)    then !y.thick = 4.0 else !y.thick = ythick
       if (n_elements(thick) eq 0)     then !p.thick = 4.0 else !p.thick = thick
       if (n_elements(charthick) eq 0) then !p.charthick = 3.0 else !p.charthick = charthick
       if (n_elements(charsize) eq 0)  then !p.charsize = 1.8 else !p.charsize = charsize
;      if (n_elements(color) eq 0)     then !p.color = fsc_color('black') else !p.color = color
    endif

    if keyword_set(keynote) then begin
       if (n_elements(xthick) eq 0)    then !x.thick = 6.0 else !x.thick = xthick
       if (n_elements(ythick) eq 0)    then !y.thick = 6.0 else !y.thick = ythick
       if (n_elements(thick) eq 0)     then !p.thick = 6.0 else !p.thick = thick
       if (n_elements(charthick) eq 0) then !p.charthick = 5.0 else !p.charthick = charthick
       if (n_elements(charsize) eq 0)  then !p.charsize = 1.8 else !p.charsize = charsize
;      if (n_elements(color) eq 0)     then !p.color = fsc_color('black') else !p.color = color
    endif

; load the plotting extras, if requested

;   if keyword_set(extras) then begin
;      plot_dir = getenv('IMPRO_DIR')+'/plot/'
;      if file_test(plot_dir,/dir) then begin
;         pushd, plot_dir
;         @myplotextras.pro
;         popd
;      endif else begin
;         splog, 'Plotting extras file '+plot_dir+'MYPLOTEXTRAS not found'
;      endelse
;   endif
    
; restore to the default IDL values

    if keyword_set(restore) then cleanplot, /silent

return
end
