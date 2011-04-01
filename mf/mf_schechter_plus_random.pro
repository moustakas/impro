function mf_schechter_plus_random, minmass, maxmass, nrand, $
  mass=mass, schechter=schechter
; jm11apr01ucsd - return a random number drawn from an double
; Schechter function

    phi = make_array(nrand,type=size(minmass,/type))
    mass = phi
    count = 0
    while (count lt nrand) do begin
       delvarx, seed1, seed2
       randmass = randomu(seed1,5000)*(maxmass-minmass)+minmass
       prob = mf_schechter_plus(randmass,schechter)

; if a second number randomly selected between 0 and 1 is less than
; that probability then keep the value, otherwise discard
       rnd = randomu(seed2,5000)
       keep = where(rnd le prob,nkeep)
       if (nkeep ne 0L) then begin
          if (count eq 0) then begin
             mass = randmass[keep]
             phi = prob[keep]
          endif else begin
             mass = [mass,randmass[keep]]
             phi = [phi,prob[keep]]
          endelse
          count = count+nkeep
       endif
    endwhile

    phi = phi[0:nrand-1]
    mass = mass[0:nrand-1]

return, phi
end
