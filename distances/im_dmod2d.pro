function im_dmod2d, dmod, sigdmod=sigdmod, sigdistance=sigdistance, A_lambda=A_lambda
; jm03may01uofa
; convert a distance modulus and error to a distance (DMOD in Mpc)
    
    ndmod = n_elements(dmod)
    if ndmod eq 0L then begin
       print, 'Syntax - distance = im_dmod2d(dmod,[sigdmod=,sigdistance=,A_lambda=])'
       return, -1L
    endif
    
    if n_elements(A_lambda) eq 0L then A_lambda = 0.0 else begin ; extinction
       if n_elements(A_lambda) gt 1L and n_elements(A_lambda) ne ndmod then begin
          print, 'A_LAMBDA and DMOD must have the same number of elements.'
          return, -1L
       endif
    endelse 

    if n_elements(sigdmod) eq 0L then sigdmod = 0.0 else begin  ; error in the distance modulus
       if n_elements(sigdmod) gt 1L and n_elements(sigdmod) ne ndmod then begin
          print, 'SIGDMOD and DMOD must have the same number of elements.'
          return, -1L
       endif
    endelse 

    dmod = dmod + A_lambda
    
    distance = 10.0^((dmod-25.0)/5.0       )     ; [Mpc]
    sigdistance = alog(10)*distance*sigdmod/5.0  ; error [Mpc]

return, distance
end
