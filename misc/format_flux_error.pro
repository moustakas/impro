pro format_flux_error, flux, ferr, newflux, newferr
; jm04apr22uofa
; jm05jul26uofa - vectorized
; given a flux and a flux error return the correct number of
; significant digits on both quantities

    nflux = n_elements(flux)
    if (nflux gt 1L) then begin
       newflux = strarr(nflux)
       newferr = strarr(nflux)
       for iflux = 0L, nflux-1L do begin
          format_flux_error, flux[iflux], ferr[iflux], newflux1, newferr1
          newflux[iflux] = strtrim(newflux1,2)
          newferr[iflux] = strtrim(newferr1,2)
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
    
return
end
    
