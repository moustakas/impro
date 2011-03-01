;+
; NAME: ARM_ZOOMPLOT
;       
; CATEGORY: plotting
;
; PURPOSE: plotting routine with zoom feature
;
; CALLING SEQUENCE: ARM_ZOOMPLOT, x, y, [coordinates=, text=, /select]
; 
; INPUTS: y - ordinate values to be plotted
;        
; OPTIONAL INPUTS: 
;   x         - abscissa values to be plotted (default=[0,1,2,3,4...])
;   text      - text message to be output for user
;   _extra    - additional parameters to be passed to DJS_PLOT
;
;   oplot_x         - abscissa values to be overplotted
;   oplot_y         - ordinate values to be overplotted
;   oplot_color     - color to be used for overplotting (string)
;   oplot_psym      - symbol to be used for overplotting
;   oplot_symsize   - symbol size to be used for overplotting   
;   oplot_thick     - thickness to be used for overplotting
;   oplot_linestyle - line style to be used for overplotting
;
;   oplotline_x         - vertical line position(s) passed to ARM_OPLOTLINE 
;   oplotline_y         - horizontal line position(s) passed to ARM_OPLOTLINE 
;   oplotline_color     - color(s) passed to ARM_OPLOTLINE
;   oplotline_thick     - thickness(es) passed to ARM_OPLOTLINE
;   oplotline_linestyle - linestyle(s) passed to ARM_OPLOTLINE
; 
;   user_xrange - default XRANGE when plotting for the first time, and
;                 when the default plot range is restored
;                 (!MOUSE.BUTTOM EQ 4)
;
; KEYWORDS:
;   silent - suppress output messages
;
; OPTIONAL OUTPUTS: coordinates - x & y coordinates of selected point
;
; PROCEDURES USED: DJS_PLOT, ARM_OPLOTLINE
;
; COMMENTS: accepts all PLOT inputs and keywords
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2003
;    improved by ARM, Oct 2003
;    middle button also return coordinates, ARM, Apr 7 2004
;    overplotting features added, ARM, May 6 2004
;    typo fixed (psym_oplot->oplot_psym), ARM, Jul 7 2004
;    typo fixed (oplotline_thick->oplot_thick), ARM, Jul 8, 2004
;    jm06jun28uofa - added USER_XRANGE optional input
; 
; Report bugs/comments to Andrew R. Marble (amarble@as.arizona.edu)
;-

pro ARM_ZOOMPLOT, xx, y, text=text, coordinates=coordinates, _extra=extra, $
                  oplot_x=oplot_x, oplot_y=oplot_y, oplot_color=oplot_color, $
                  oplot_psym=oplot_psym, oplot_linestyle=oplot_linestyle, $
                  oplot_thick=oplot_thick, oplot_symsize=oplot_symsize, $
                  oplotline_x=oplotline_x, oplotline_color=oplotline_color, $
                  oplotline_y=oplotline_y, oplotline_thick=oplotline_thick, $
                  oplotline_linestyle=oplotline_linestyle, single=single, $
                  silent=silent, user_xrange=user_xrange

; define defaults

    if N_ELEMENTS(y) eq 0L then begin
       y = xx
       x = INDGEN(N_ELEMENTS(y))
    endif else x = xx
    if N_ELEMENTS(y) ne N_ELEMENTS(x) then $
      MESSAGE, 'X and Y have incompatible dimensions.'

    if N_ELEMENTS(oplot_y) ne 0L then begin
       if N_ELEMENTS(oplot_x) eq 0L then oplot_x=INDGEN(N_ELEMENTS(oplot_y))
       if N_ELEMENTS(oplot_color) eq 0L then oplot_color=''
       if N_ELEMENTS(oplot_psym) eq 0L then oplot_psym=0
       if N_ELEMENTS(oplot_symsize) eq 0L then oplot_symsize=1
       if N_ELEMENTS(oplot_linestyle) eq 0L then oplot_linestyle=0
       if N_ELEMENTS(oplot_thick) eq 0L then oplot_thick=1
    endif
    if N_ELEMENTS(oplot_y) ne N_ELEMENTS(oplot_x) then $
      MESSAGE, 'OPLOT_X and OPLOT_Y have incompatible dimensions.'

    if N_ELEMENTS(oplotline_x) + N_ELEMENTS(oplotline_y) ne 0L then begin
       if N_ELEMENTS(oplotline_color) eq 0L then oplotline_color=''
       if N_ELEMENTS(oplotline_thick) eq 0L then oplotline_thick=1
       if N_ELEMENTS(oplotline_linestyle) eq 0L then oplotline_linestyle=1
    endif

; copy original abscissa values

    xcopy = x

; plot x versus y

    DJS_PLOT, x, y, xrange=user_xrange, _extra=extra
    if N_ELEMENTS(oplotline_x) ne 0L then ARM_OPLOTLINE, oplotline_x, $
      color=oplotline_color, thick=oplotline_thick, linestyle=oplotline_linestyle
    if N_ELEMENTS(oplotline_y) ne 0L then ARM_OPLOTLINE, oplotline_y, $
      color=oplotline_color, thick=oplotline_thick, linestyle=oplotline_linestyle, $
      /horizontal
    if N_ELEMENTS(oplot_y) ne 0L then DJS_OPLOT, oplot_x, oplot_y, $
      color=oplot_color, psym=oplot_psym, linestyle=oplot_linestyle, $
      thick=oplot_thick, symsize=oplot_symsize

; output display message

    if not KEYWORD_SET(silent) then begin
       if not KEYWORD_SET(text) then $
         text = 'Use mouse to zoom in/out on plot...'
       PRINT, text
       PRINT
       PRINT, '   left button: select opposite corners of zoom region'
       PRINT, ' middle button: print cursor coordinates (double-click to exit)'
       PRINT, '  right button: restore original plotting range'
    endif

    !mouse.button = 0           ; initialize mouse selection

    done = 0L
    oldx = 9d9 & oldy = -9d9

    while not done do begin

       CURSOR, xcoord1, ycoord1, 3, /data ; select mouse position

       if !mouse.button eq 4 then begin

; restore original plot range

          ERASE
          DJS_PLOT, x, y, xrange=user_xrange, /noerase, _extra=extra
          if N_ELEMENTS(oplotline_x) ne 0L then ARM_OPLOTLINE, oplotline_x, $
            color=oplotline_color, thick=oplotline_thick, $
            linestyle=oplotline_linestyle
          if N_ELEMENTS(oplotline_y) ne 0L then ARM_OPLOTLINE, oplotline_y, $
            color=oplotline_color, thick=oplotline_thick, $
            linestyle=oplotline_linestyle, /horizontal
          if N_ELEMENTS(oplot_y) ne 0L then DJS_OPLOT, oplot_x, oplot_y, $
            color=oplot_color, psym=oplot_psym, linestyle=oplot_linestyle, $
            thick=oplot_thick, symsize=oplot_symsize
       

       endif else if !mouse.button eq 1 then begin

; select opposite corner of zoom region

          !mouse.button = 0 ; reset mouse selection

          while !mouse.button ne 1 do begin

             CURSOR, xcoord2, ycoord2, 3, /data 
             if !mouse.button ne 1 then message, /continue, $
               'use LEFT mouse button to select zoom region'

          endwhile

; zoom in on desired region

          ERASE
          DJS_PLOT, x, y, xrange=MINMAX([xcoord1,xcoord2]), xstyle=1, $
            yrange=MINMAX([ycoord1,ycoord2]), ystyle=1, /noerase, _extra=extra
          if N_ELEMENTS(oplotline_x) ne 0L then ARM_OPLOTLINE, oplotline_x, $
            color=oplotline_color, thick=oplotline_thick, linestyle=oplotline_linestyle
          if N_ELEMENTS(oplotline_y) ne 0L then ARM_OPLOTLINE, oplotline_y, $
            color=oplotline_color, thick=oplotline_thick, /horizontal, $
            linestyle=oplotline_linestyle
          if N_ELEMENTS(oplot_y) ne 0L then DJS_OPLOT, oplot_x, oplot_y, $
            color=oplot_color, psym=oplot_psym, linestyle=oplot_linestyle, $
            thick=oplot_thick, symsize=oplot_symsize

       endif else begin
          
          if keyword_set(single) or (xcoord1 eq oldx and ycoord1 eq oldy) then begin
;             coordinates = [xcoord1, ycoord1]
             done = 1L 
          endif else begin
             print, xcoord1, ycoord1
             oldx = xcoord1
             oldy = ycoord1
          endelse
          
       endelse
       
       coordinates = [xcoord1, ycoord1]

    endwhile
    
    x = xcopy ; restore unaltered values

    return
    
 end
 
