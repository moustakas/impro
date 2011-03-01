;+
; NAME: arm_winnow
;       
; CATEGORY: plotting
;
; PURPOSE: select/unselect data points on a plot
;
; CALLING SEQUENCE: result = arm_winnow(x, y, [mask=, /quiet, /single,
;                                              extra=] 
;
; INPUTS:
;   x - array of plotted abscissa values
;   y - array of plotted ordinate values
;       
; OPTIONAL INPUTS:
;   mask  - array of selected values (0=unselected, 1=selected; all
;           selected by default) 
;   extra - paramaters passed to DJS_PLOT (for matching existing plot)
;           (eg, if the original plot had yellow triangles, then
;           passing PSYM=5 and COLOR='yellow' would result in points
;           which are unselected and then reselected being restored to
;           their original appearance)
;
; KEYWORDS:
;   single - return after a single selection
;   quiet  - suppress instructions
;
; OUTPUTS:
;   returns updated mask array (0=unselected, 1=selected)
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 April 9
;-

function arm_winnow, x, y, mask=mask, quiet=quiet, psym=psym, color=color, single=single

    n = N_ELEMENTS(x) ; number of data points

; defaults and error checking

    if N_ELEMENTS(mask) eq 0L then mask = LONARR(n) + 1L
    if N_ELEMENTS(mask) ne n or N_ELEMENTS(y) ne n then $
      MESSAGE, 'X, Y and MASK have incompatible dimensions!'
    if N_ELEMENTS(psym) eq 0L or N_ELEMENTS(psym) gt 1L then psym = 2

; print instructions if desired

    if not KEYWORD_SET(quiet) then begin
       PRINT
       PRINT, 'LEFT   mouse button UNSELECTS nearest selected point.'
       PRINT, 'MIDDLE mouse button SELECTS nearest unselected point.'
       PRINT, 'RIGHT  mouse button exits.'
       PRINT
    endif

; convert data coordinates to device coordinates

    coords = convert_coord(x, y, /data, /to_device)
    xx = REFORM(coords[0,*])
    yy = REFORM(coords[1,*])

    !mouse.button = 0L ; initialize mouse input

; loop until right mouse button is pressed
    
    while !mouse.button ne 4L do begin

       CURSOR, xmouse, ymouse, 3, /device ; read mouse input

       if !mouse.button ne 4L then begin

          distance = SQRT((xx-xmouse)^2 + (yy-ymouse)^2)
          selected = WHERE(mask eq 1L, complement=unselected, nselected)
          
; select or unselect nearest data point

          if !mouse.button eq 1L then begin
             if nselected gt 0L then begin
                wh = (WHERE(distance[selected] eq MIN(distance[selected])))[0]
                mask[selected[wh]] = 0L
                DJS_OPLOT, x[selected[wh]]*[1,1], y[selected[wh]]*[1,1], psym=7, thick=2, color='red'
             endif else MESSAGE, 'All points are already unselected!', /continue
          endif
          if !mouse.button eq 2L then begin
             if nselected lt n then begin
                wh = (WHERE(distance[unselected] eq MIN(distance[unselected])))[0]
                mask[unselected[wh]] = 1L
                OPLOT, x[unselected[wh]]*[1,1], y[unselected[wh]]*[1,1], psym=7, thick=2, color=!p.background
                DJS_OPLOT, x[unselected[wh]]*[1,1], y[unselected[wh]]*[1,1], extra=extra
             endif else MESSAGE, 'All points are already selected!', /continue
          endif
          
       endif
       
       if KEYWORD_SET(single) then !mouse.button = 4L
       
    endwhile
    
    return, mask
    
 end
