;+
; NAME:
;   FORMAT_FLUX_ERROR
;
; PURPOSE:
;   Silly routine to truncate a vector to the correct number of
;   significant digits, given the uncertainty on that measurement. 
;
; INPUTS: 
;   flux - input vector (of any kind, not necessarily flux)
;   ferr - error on FLUX
;
; OUTPUTS: 
;   newflux - string output with correct number of significant digits 
;   newferr - corresponding string error on NEWFLUX
;
; COMMENTS:
;   This should be generalized and expanded, but it's right and
;   it works. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Apr 22, U of A
;   jm05jul26uofa - vectorized
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

pro format_flux_error, flux, ferr, newflux, newferr

    nflux = n_elements(flux)
    if (nflux gt 1L) then begin
       newflux = strarr(nflux)
       newferr = strarr(nflux)
       for iflux = 0L, nflux-1 do begin
          format_flux_error, flux[iflux], ferr[iflux], newflux1, newferr1
          newflux[iflux] = newflux1
          newferr[iflux] = newferr1
       endfor
       return
    endif
    
    if (flux ge 0.01) and (flux lt 0.1) then begin
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,3),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
    endif
    
    if (flux ge 0.1) and (flux lt 1.0) then begin
       if (ferr ge 0.001) and (ferr lt 0.01) then begin
          newflux = string(fix_digits(flux,4),format='(F12.4)')
          newferr = string(fix_digits(ferr,4),format='(F12.4)')
       endif
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,3),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,2),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
    endif
    
    if (flux ge 1.0) and (flux lt 10.0) then begin
       if (ferr ge 0.001) and (ferr lt 0.01) then begin
          newflux = string(fix_digits(flux,4),format='(F12.4)')
          newferr = string(fix_digits(ferr,3),format='(F12.4)')
       endif
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,4),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,3),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,2),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
    endif
    
    if (flux ge 10.0) and (flux lt 100.0) then begin
       if (ferr ge 0.001) and (ferr lt 0.01) then begin
          newflux = string(fix_digits(flux,6),format='(F12.4)')
          newferr = string(fix_digits(ferr,4),format='(F12.4)')
       endif
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,5),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,4),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,3),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
       if (ferr ge 10.0) and (ferr lt 100.0) then begin
          newflux = string(fix_digits(flux,2),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
    endif
    
    if (flux ge 100.0) and (flux lt 1000.0) then begin
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,6),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,5),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,4),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
       if (ferr ge 10.0) and (ferr lt 100.0) then begin
          newflux = string(fix_digits(flux,3),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 100.0) and (ferr lt 1000.0) then begin
          newflux = string(fix_digits(flux,2),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
    endif

    if (flux ge 1000.0) and (flux lt 10000.0) then begin
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,7),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,6),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,5),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
       if (ferr ge 10.0) and (ferr lt 100.0) then begin
          newflux = string(fix_digits(flux,4),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 100.0) and (ferr lt 1000.0) then begin
          newflux = string(fix_digits(flux,3),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 1000.0) and (ferr lt 10000.0) then begin
          newflux = string(fix_digits(flux,2),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
    endif

    if (flux ge 10000.0) and (flux lt 100000.0) then begin
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,8),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,7),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,6),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
       if (ferr ge 10.0) and (ferr lt 100.0) then begin
          newflux = string(fix_digits(flux,5),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 100.0) and (ferr lt 1000.0) then begin
          newflux = string(fix_digits(flux,4),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 1000.0) and (ferr lt 10000.0) then begin
          newflux = string(fix_digits(flux,3),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 10000.0) and (ferr lt 100000.0) then begin
          newflux = string(fix_digits(flux,2),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
    endif

    if (flux ge 100000.0) and (flux lt 1000000.0) then begin
       if (ferr ge 0.01) and (ferr lt 0.1) then begin
          newflux = string(fix_digits(flux,9),format='(F12.3)')
          newferr = string(fix_digits(ferr,3),format='(F12.3)')
       endif
       if (ferr ge 0.1) and (ferr lt 1.0) then begin
          newflux = string(fix_digits(flux,8),format='(F12.2)')
          newferr = string(fix_digits(ferr,2),format='(F12.2)')
       endif
       if (ferr ge 1.0) and (ferr lt 10.0) then begin
          newflux = string(fix_digits(flux,7),format='(F12.1)')
          newferr = string(fix_digits(ferr,2),format='(F12.1)')
       endif
       if (ferr ge 10.0) and (ferr lt 100.0) then begin
          newflux = string(fix_digits(flux,6),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 100.0) and (ferr lt 1000.0) then begin
          newflux = string(fix_digits(flux,5),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 1000.0) and (ferr lt 10000.0) then begin
          newflux = string(fix_digits(flux,4),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 10000.0) and (ferr lt 100000.0) then begin
          newflux = string(fix_digits(flux,3),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
       if (ferr ge 100000.0) and (ferr lt 1000000.0) then begin
          newflux = string(fix_digits(flux,2),format='(I12)')
          newferr = string(fix_digits(ferr,2),format='(I12)')
       endif
    endif

    if (n_elements(newflux) eq 0L) then message, 'Code me up!'

    newflux = strtrim(newflux,2)
    newferr = strtrim(newferr,2)
    
return
end
    
