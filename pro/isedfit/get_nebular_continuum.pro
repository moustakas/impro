;+
; NAME:
;   GET_NEBULAR_CONTINUUM()
;
; PURPOSE:
;   Given the number of Lyman-continnum photons, compute the nebular
;   continuum spectrum.
;
; INPUTS:
;   nlyc - number of Lyman-continuum photons [sec^-1]
;
; OPTIONAL INPUTS:
;   wave - input/output rest-frame wavelength vector [Angstrom]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   flam_neb - nebular continuum spectrum in units of erg/s/A [NSPEC] 
;
; COMMENTS:
;   See Koleva et al. (2009, Section 3.3) for many details. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 10, Siena
;
; Copyright (C) 2013, John Moustakas
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

; references: Leitherer & Heckman 95; Koleva+09; 

Function get_nebular_continuum, nlyc, wave=wave, temp=temp

    if n_elements(nlyc) eq 0 then begin
       doc_library, 'get_nebular_continuum'
       return, -1
    endif

; build the wavelength vector and the output spectrum 
    if n_elements(wave) eq 0 then wave = range(912.0,1E10,500) ; [Angstrom]
    nwave = n_elements(wave)
    flam_neb = wave*0.0

; total recombination coefficient alpha2(T)
    alpha2 = 2.575D-13 ; Aller+84 [cm^3 s^-1]
;   alpha2 = 2.616D-13 ; Pegase.2

; assumed densities of He+ and He++ relative to H+, taken from
; Koleva+09 
    f_HeI = 0.0897    ; Pegase uses 0.095
    f_HeII = 1.667D-4 ; Pegase assumes zero

    light = im_light(/ang) ; [Angstrom/s]

; eventually replace this bit of code with my own computations of the
; total emission coefficient (see for example Koleva+09), but for now
; just use Pegase's calculations
    peg = pegase_hii()
    
         fneb(i)=1.d-40*(g1+g2+HeI*g3+HeII*g4)*c/alphaTe        /lambda(i)**2
  
    
; total emission coefficient
    gamma_tot = gamma_b + gamma_2p + gamma_HI + $
      f_HeI*gamma_HeI + f_HeII*gamma_HeII

; see equation (1) in Koleva+09    
    flam_neb = gamma_tot*nlyc*light/wave^2/alpha2
    
return, flam_neb
end
