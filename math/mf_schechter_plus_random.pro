;+
; NAME:
;   MF_SCHECHTER_PLUS_RANDOM()
;
; PURPOSE:
;   Draw random numbers from a (specified) double Schechter function.  
;
; INPUTS: 
;   schechter - MF_SCHECHTER_PLUS() style structure:
;     .PHISTAR - number density at the 'knee' of the MF
;     .LOGMSTAR - log-base-10 stellar mass at the 'knee' of the MF
;     .ALPHA - low-mass slope
;     .PHIPLUS - second normalization at low mass
;     .ALPHAPLUS - second slope at low mass
;
; OPTIONAL INPUTS: 
;   minmass - minimum stellar mass to simulate, log-base-10 in units
;     of Msun (default 8)
;   maxmass - maximum stellar mass to simulate (default 12)
;   nrand - number of random numbers to return (default 1000)
;
; KEYWORD PARAMETERS: 
;   sortit - sort the output by mass
;
; OUTPUTS: 
;   phi - Phi(log M) mass function values corresponding to MASS
;     [NRAND] 
;
; OPTIONAL OUTPUTS:
;   mass - random stellar masses [NRAND]
;
; COMMENTS:
;   Note that Phi(log M) = ln(10)*M*Phi(M)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Apr 01, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function mf_schechter_plus_random, schechter, minmass=minmass, $
  maxmass=maxmass, nrand=nrand, mass=mass, sortit=sortit

    if (n_elements(schechter) eq 0) then begin
       doc_library, 'mf_schechter_random'
       return, -1
    endif

    if (tag_exist(schechter,'phistar') eq 0) or $
      (tag_exist(schechter,'logmstar') eq 0) or $
      (tag_exist(schechter,'alpha') eq 0) or $
      (tag_exist(schechter,'phiplus') eq 0) or $
      (tag_exist(schechter,'alphaplus') eq 0) then begin
       splog, 'Improper double SCHECHTER structure!'
       return, -1
    endif

    if (n_elements(minmass) eq 0) then minmass = 8.0  ; [Msun]
    if (n_elements(maxmass) eq 0) then maxmass = 12.0 ; [Msun]
    if (n_elements(nrand) eq 0) then nrand = 1000
    
    phi = lonarr(nrand)
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

; truncate to the desired number of random numbers
    phi = phi[0:nrand-1]
    mass = mass[0:nrand-1]

    if keyword_set(sortit) then begin
       srt = sort(mass)
       phi = phi[srt]
       mass = mass[srt]
    endif
    
return, phi
end
