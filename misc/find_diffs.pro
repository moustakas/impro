function find_diffs, x, y=y
; jm01jul23uofa
; return all differences of a vector x, or between two vectors x and y 

    nx = n_elements(x)

    if keyword_set(y) then begin

       ny = n_elements(y)
       ndiffs = nx*ny
       diffs = fltarr(ndiffs)

       count = 0L
       for i = 0L, nx-1L do begin

          for j = 0L, ny-1L do begin
             diffs[count] = abs(x[i]-y[j])
             count = count+1L
          endfor
       
       endfor

    endif else begin
    
       ndiffs = nx*(nx-1L)/2L
       diffs = fltarr(ndiffs)

       count = 0L
       for i = 0L, nx-2L do begin

          for j = i+1L, nx-1L do begin
             diffs[count] = abs(x[i]-x[j])
             count = count+1L
          endfor
       
       endfor

    endelse

return, diffs
end
