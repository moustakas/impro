;+
; NAME:
;       COMPUTE_LINELUMS()
;
; PURPOSE:
;       Compute emission-line luminosities and other distant-dependent
;       quantities.  
;
; INPUTS:
;       linefit - data structure from ISPECLINEFIT()
;
; OPTIONAL INPUTS:
;       ancillary  - ancillary data (specifically, a data structure
;                    containing the distance)
;       disttag    - name of the distance tag name (default DISTANCE)  
;       disterrtag - name of the distance error tag name (default
;                    DISTANCE_ERR)  
;       photerrtag - name of the absolute photometric error tag name
;                    (the default is to assume that this tag does not
;                    exist)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       linelums - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 11, U of A
;       jm04jan12uofa - add the R-band magnitude
;       jm04feb06uofa - major updates, remove magnitude calculations, 
;                       compute SFR
;       jm04mar11uofa - check for the distance in LINEFIT
;       jm04jul22uofa - added _EXTRA keyword
;       jm05jan07uofa - removed the SFR computations
;                       to a new routine, COMPUTE_SFRS() 
;       jm06apr11uofa - minor cleanup
;       jm08feb05nyu  - cleaned up 
;
; Copyright (C) 2003-2006, 2008, John Moustakas
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

function compute_linelums, linefit, disttag=disttag, disterrtag=disterrtag, $
  photerrtag=photerrtag, ancillary=ancillary, select_lines=select_lines, $
  _extra=extra

    lsun = 3.826D33    ; [erg/s]
    Mpc2cm = 3.086D24  ; [cm/Mpc]
    lsunerg = 3.826D33 ; erg/s

    nspec = n_elements(linefit)
    if (nspec eq 0L) then begin
       doc_library, 'compute_linelums'
       return, -1L
    endif

    nselect = n_elements(select_lines)
    if (nselect eq 0L) then begin
       splog, 'SELECT_LINES input required.'
       return, -1L
    endif

    nancillary = n_elements(ancillary)
    if (nancillary eq 0L) then ancillary = linefit

; initialize the output data structure    
    
    linelums = create_struct(select_lines[0]+'_lum',[-999.0,-999.0])
    for k = 1L, nselect-1L do $
      linelums = create_struct(linelums,select_lines[k]+'_lum',[-999.0,-999.0])    
    linelums = replicate(struct_addtags({linelumname: select_lines},linelums),nspec)
    
; identify the appropriate tag names, etc.
    
    if (n_elements(disttag) eq 0L) then disttag = 'DISTANCE'
    if (n_elements(disterrtag) eq 0L) then disterrtag = 'DISTANCE_ERR'

    if (tag_exist(ancillary,disttag) eq 0) then begin
       splog, 'No DISTTAG field in LINEFIT.'
       return, linelums
    endif else begin
       match = where(disttag eq tag_names(ancillary),nmatch)
       if (nmatch eq 1L) then distance = ancillary.(match) else begin
          splog, 'Problem finding DISTTAG field in ANCILLARY.'
          return, linelums
       endelse
    endelse
    
    if (tag_exist(ancillary,disterrtag) eq 0) then begin
;      splog, 'No DISTERRTAG field in LINEFIT!'
       distance_err = distance*0.0 ; assumed no error
    endif else begin
       match = where(disterrtag eq tag_names(ancillary),nmatch)
       if (nmatch eq 1L) then distance_err = ancillary.(match) else begin
          splog, 'Problem finding DISTERRTAG field in ANCILLARY.'
          return, linelums
       endelse
    endelse

    if (n_elements(photerrtag) ne 0L) then begin
       if (tag_exist(ancillary,photerrtag) eq 0) then begin
;         splog, 'No PHOTERRTAG field in ANCILLARY.'
          photerror = replicate(0.0,nspec) ; [%]
       endif else begin
          match = where(photerrtag eq tag_names(ancillary),nmatch)
          if (nmatch eq 1L) then photerror = ancillary.(match) else begin
             splog, 'Problem finding PHOTERRTAG field in ANCILLARY.'
             photerror = distance*0.0 ; [%]
          endelse
       endelse
    endif else photerror = distance*0.0
    
    for j = 0L, nselect-1L do begin

       mini = struct_trimtags(linefit,select=select_lines[j])

       array = mini.(0)
       flux = reform(array[0,*])
       flux_err = reform(array[1,*])

       well = where(flux_err gt 0.0,nwell)
       if (nwell ne 0L) then flux_err[well] = sqrt( flux_err[well]^2.0 + (photerror*flux[well]/100.0)^2.0 )
       
       good = where((distance gt -900.0) and (flux/flux_err gt 0.0),ngood)
       if (ngood ne 0L) then begin

          dist = distance[good]*Mpc2cm
          nodisterr = where(distance_err[good] lt -900.0,nnodisterr)
          if nnodisterr ne 0L then distance_err[good] = 0.0
          dist_err = distance_err[good]*Mpc2cm

          flux = flux[good]
          flux_err = sqrt( (flux_err[good])^2.0 + (photerror[good]*flux/100.0)^2.0 )
          
          lum = flux*4.0*!dpi*dist*dist                                                ; [erg/s]
          lum_err = 4.0*!dpi*sqrt((2*flux*dist*dist_err)^2.0+(dist*dist*flux_err)^2.0) ; [erg/s]

; take the logarithm and fill the structure       

          lumdata = transpose([ [alog10(lum/lsun)], [lum_err/lum/alog(10.0)] ]) ; luminosity and error
          linelums[good].(j+1) = lumdata 

;         linelums[good].(2*j+1) = lum_err / lum / alog(10.0) ; luminosity error
;         linelums[good].(2*j) = alog10(lum/lsun)             ; luminosity

       endif
       
    endfor

return, linelums
end    
