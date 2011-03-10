function asinh_random, x1, x2, num, soft=soft
; jm11mar08ucsd - return a random number variable drawn from an
; ASINH() distribution (see RANGE); deals with multidimensional arrays
; (although it's pretty slow)

    if (n_elements(soft) gt 0) then ss = soft else ss = 0.3D

    out = make_array(num,type=size(x1,/type))
    count = 0
    while (count lt cmproduct(num)) do begin
       delvarx, seed1, seed2
       xx = randomu(seed1,1)
       sh = sinh(xx[0]/ss)
       prob = x1 + (x2-x1)*sh/(sinh(1.0/ss))

; if a second number randomly selected between 0 and 1 is less than
; that probability then keep the value, otherwise discard
       rnd = randomu(seed2,1)
       if (rnd[0] le prob) then begin
          out[count] = prob
          count++
       endif
    endwhile

return, out
end
