function maggies2mag, maggies, ivarmaggies=ivarmaggies, magerr=magerr
; jm09may24nyu - convert maggies to AB magnitude; also optionally 
;   converts inverse variance maggies to (AB) magnitude error
    mag = maggies*0.0
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
    
