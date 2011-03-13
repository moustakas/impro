function mag2maggies, mag, magerr=magerr, vega2ab=vega2ab, $
  ivarmaggies=ivarmaggies
; jm09may24nyu - convert AB (or Vega) magnitudes to maggies and, given
;   the magnitude error, to inverse variance maggies

    ndim = size(mag,/n_dim)
    dim = size(mag,/dim)
    nband = dim[0]
    if (ndim eq 1) then nobj = 1 else nobj = dim[1]
    
    if (n_elements(vega2ab) eq 0) then vega2ab = fltarr(nband)
    maggies = mag*0.0
    ivarmaggies = mag*0.0
    if (n_elements(magerr) eq 0) then doivar = 0 else doivar = 1
    for ii = 0, nband-1 do begin
       good = where((mag[ii,*] gt 0.0) and (mag[ii,*] lt 90.0),ngood)
       if (ngood ne 0L) then begin
          mag1 = mag[ii,good] + vega2ab[ii]
          maggies[ii,good] = 10.0^(-0.4*mag1)
       endif
       if doivar then begin
          good = where((mag[ii,*] gt 0.0) and (mag[ii,*] lt 90.0) and $
            (magerr[ii,*] gt 0.0),ngood)
          if (ngood ne 0L) then begin
             ivarmaggies[ii,good] = 1.0/(0.4*alog(10.0)*(maggies[ii,good]*magerr[ii,good]))^2
          endif
       endif
    endfor

return, maggies
end
    
