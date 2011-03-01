function mag2maggies, mag, magerr=magerr, vega2ab=vega2ab, $
  ivarmaggies=ivarmaggies
; jm09may24nyu - convert AB (or Vega) magnitudes to maggies and, given
;   the magnitude error, to inverse variance maggies
    dim = size(mag,/dim)
    nband = dim[0]
    nobj = dim[1]
    
;   if (n_elements(vega2ab) eq 0) then vega2ab = fltarr(nband)
;   maggies = mag*0.0
;   for ii = 0, nband-1 do begin
;      good = where((mag[ii,*] gt 0.0) and (mag[ii,*] lt 90.0)
;      mag1 = mag[ii,*] + vega2ab[ii]
;      
;      maggies[ii
;   endfor

    good = where(maggies gt 0.0,ngood)
    if (ngood ne 0L) then mag[good] = -2.5*alog10(maggies[good])
    if arg_present(magerr) and (n_elements(ivarmaggies) ne 0) then begin
       magerr = ivarmaggies*0.0
       good = where((maggies gt 0.0) and (ivarmaggies gt 0.0),ngood)
       if (ngood ne 0L) then magerr[good] = 1.0/(0.4*sqrt(ivarmaggies[good])*$
         maggies[good]*alog(10.0))
    endif
return, mag
end
    
