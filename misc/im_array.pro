function im_array, mmin, mmax, dm, double=double, $
  narr=narr, log=log
; jm09mar02nyu - build a simple linear or logarithmic array
; jm09oct08ucsd - added LOG keyword

    if (n_elements(mmin) eq 0L) then mmin = 0.0
    if (n_elements(mmax) eq 0L) then mmax = 1.0
    if (n_elements(dm) eq 0L) then dm = 0.1

    if keyword_set(log) then begin
       logarr = dindgen((alog10(mmax)-alog10(mmin))/dm+1)*dm+alog10(mmin)
       arr = 10.0^logarr
    endif else arr = dindgen((mmax-mmin)/dm+1)*dm+mmin
    narr = n_elements(arr)

    if (keyword_set(double) eq 0) then arr = float(arr)
    
return, arr
end
