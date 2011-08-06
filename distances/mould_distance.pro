;+
; NAME:
;       MOULD_DISTANCE()
;
; PURPOSE:
;       Compute the Hubble-flow velocity and luminosity distance of an
;       object corrected for infall onto Virgo, the Great Attractor,
;       and the Shapley concentration using the methodology described
;       in Mould et al. (2000).
;
; CALLING SEQUENCE:
;       mould = mould_distance(ra,dec,cz,object=,H0=,$
;          omega0=,omega_lambda=,/proper)
;
; INPUTS:
;       ra  - right ascension (J2000) [HMS]
;       dec - declination (J2000) [DMS]
;       cz  - heliocentric radial velocity [km/s]
;
; OPTIONAL INPUTS:
;       object       - optional object name
;       H0           - Hubble constant (default 70 km/s/Mpc)
;       omega0       - matter density (default 0.3)
;       omega_lambda - cosmological constant density (default 0.7)
;
; KEYWORD PARAMETERS:
;       proper - return the proper distance (rather than the
;                luminosity distance)
;
; OUTPUTS:
;       mould - output data structure with the following fields: 
;          object   - object name
;          ra       - input declination
;          dec      - input right ascension
;          z        - input heliocentric redshift
;          cz       - input heliocentric velocity distance
;          v_LG     - Local Group velocity [km/s]
;          v_infall - infall velocity [km/s] (three components)
;          v_cosmic - Hubble-flow velocity [km/s]
;          distance - luminosity distance [Mpc]
;          flag     - flag [0=good, 1=Virgo, 2=GA, 3=Shapley]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       See Mould et al. 2000, ApJ, 529, 786 for the formalism. 
;
; PROCEDURES USED:
;       GLACTC, DJS_DIFF_ANGLE(), IM_HMS2DEC()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 October 13, U of A - written 
;       jm04apr12uofa - minor bug fixes; properly zero the velocity to
;                       the attractor center
;       jm04may19uofa - updated and cleaned up documentation 
;
; Copyright (C) 2003-2004, John Moustakas
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

function mould_distance, ra, dec, cz, object=object, H0=H0, $
  omega0=omega0, omega_lambda=omega_lambda, proper=proper

    light = 2.99792458D5 ; speed of light [km/s]
    
    nra = n_elements(ra)
    ndec = n_elements(dec)
    ncz = n_elements(cz)
    nobject = n_elements(object)
    
    if (nra eq 0L) or (ndec eq 0L) or (ncz eq 0L) then begin
       print, 'Syntax - mould = mould_distance(ra,dec,cz,object=,$'
       print, '   H0=,omega0=,omega_lambda=,/proper)'
       return, -1L
    endif

    if (nra ne ndec) then begin
       print, 'RA and DEC have different numbers of elements.'
       return, -1L
    endif

    if (nra ne ncz) then begin
       print, 'RA and CZ have different numbers of elements.'
       return, -1L
    endif

    if nobject eq 0L then begin
       object = replicate('',nra)
       nobject = nra
    endif

    if (nra ne nobject) then begin
       print, 'RA and OBJECT have different numbers of elements.'
       return, -1L
    endif

    if n_elements(H0) eq 0L then H0 = 70.0 ; [km/s/Mpc]
    if n_elements(omega0) eq 0L then omega0 = 0.3
    if n_elements(omega_lambda) eq 0L then omega_lambda = 0.7

    q0 = omega0/2.0 - omega_lambda

    d_Virgo = 15.3 ; Virgo luminosity distance [Mpc]
    
; Table A1; references:  Aaronson et al (1982), Han (1992), Faber &
; Burstein (1989), Shaya, Tully, & Pierce (1992), and Huchra (1995)
    
    model = {$
      cluster:    '',   $
      ra:         '',   $
      dec:        '',   $
      ra_deg:     0.0D, $ ; galaxy right ascension
      dec_deg:    0.0D, $ ; galaxy declination
      v_helio:    0.0D, $ ; observed (mean) heliocentric velocity
      v_LG:       0.0D, $ ; velocity corrected to the centroid of the LG
      v_fid:      0.0D, $ ; adopted model infall velocity at the position of the LG
      radius:     0.0D, $ ; cluster radius [degrees]
      v_range_lo: 0.0D, $
      v_range_hi: 0.0D}
    model = replicate(model,3)

    model.cluster = ['Virgo','GA','Shapley']
    model.ra      = ['12:28:19','13:20:00','13:30:00']
    model.dec     = ['+12:40:00','-44:00:00','-31:00:00']
    model.ra_deg  = 15.0D*im_hms2dec(model.ra)
    model.dec_deg = im_hms2dec(model.dec)
    model.v_helio = [1035,4600,13800]     ; [km/s]
    model.v_LG    = [957,4380,13600]      ; [km/s]
    model.v_fid   = [200,400,85]          ; [km/s]
    model.radius  = [10,10,12]            ; [degrees]
    model.v_range_lo = [600,2600,10000]   ; [km/s]
    model.v_range_hi = [2300,6600,16000]  ; [km/s]
    
    N = 3       ; number of concentrations (Virgo, GA, Shapley)
    gamma = 2.0 ; mass concentration density profile

; initialize the output data structure

    mould = {$
      object:                   '', $ ; object name
      ra:                       '', $ ; input declination
      dec:                      '', $ ; input right ascension
      z:                      0.0D, $ ; input heliocentric redshift
      cz:                     0.0D, $ ; input heliocentric velocity distance
      v_LG:                  -999D, $ ; Local Group velocity [km/s]
      v_infall: replicate(-999D,N), $ ; infall velocity [km/s] (three components)
      v_cosmic:              -999D, $ ; Hubble-flow velocity [km/s]
      distance:              -999D, $ ; luminosity distance [Mpc]
      flag:                      0L}  ; flag [0=good, 1=Virgo, 2=GA, 3=Shapley]
    mould = replicate(mould,nobject)

    mould.object = object
    mould.ra = ra
    mould.dec = dec
    mould.cz = cz
    
    mould.z = cz/light

; loop on each object

    for iobj = 0L, nobject-1L do begin

; begin with false flags
       
       virgo = 0
       GA = 0
       shapley = 0
       
; convert to Galactic coordinates

       ra_deg = 15.0D*im_hms2dec(ra[iobj]) ; [degrees]
       dec_deg = im_hms2dec(dec[iobj])     ; [degrees]
       
       glactc, ra_deg, dec_deg, 2000.0, gl, gb, 1, /degree
       gl = gl*!dtor
       gb = gb*!dtor

; first convert the galaxy velocity to the Local Group frame using the
; prescription in Yahil, Tammann, & Sandage (1977) (equation A1); we
; assume that to first order that the apparent radial velocity of an
; object in the Local Group frame represents its distance

       mould[iobj].v_LG = cz[iobj] - 79.0*cos(gl)*cos(gb) + 296.0*sin(gl)*cos(gb) - 36.0*sin(gb)
       
; compute the angular and velocity object-attractor distances
; (equation 2)
       
       theta = djs_diff_angle(replicate(ra_deg,N),replicate(dec_deg,N),$
         model.ra_deg,model.dec_deg,units='degrees')*!dtor

       r0a = sqrt(mould[iobj].v_LG^2.0 + model.v_LG^2.0 - 2*mould[iobj].v_LG*model.v_LG*cos(theta))

; check to see if this object is in the direction of Virgo, the Great
; Attractor, or the Shapley concentration and assign the proper flag
; and distance; ; otherwise compute the infall velocity correction and
; the corresponding cosmic velocity (equations 1 & A2)

       if (theta[0]*!radeg lt model[0].radius) and ((model[0].v_LG-r0a[0] gt model[0].v_range_lo) or $
         (model[0].v_LG+r0a[0] lt model[0].v_range_hi)) then virgo = 1

       if (theta[1]*!radeg lt model[1].radius) and ((model[1].v_LG-r0a[1] gt model[1].v_range_lo) or $
         (model[1].v_LG+r0a[1] lt model[1].v_range_hi)) then GA = 1

       if (theta[2]*!radeg lt model[2].radius) and ((model[2].v_LG-r0a[2] gt model[2].v_range_lo) or $
         (model[2].v_LG+r0a[2] lt model[2].v_range_hi)) then shapley = 1

       if virgo or GA or shapley then begin

          mould[iobj].v_infall = 0.0

          if virgo then begin
             mould[iobj].v_cosmic = model[0].v_LG
             mould[iobj].flag = 1L
          endif
          if GA then begin
             mould[iobj].v_cosmic = model[1].v_LG
             mould[iobj].flag = 2L
          endif
          if shapley then begin
             mould[iobj].v_cosmic = model[2].v_LG
             mould[iobj].flag = 3L
          endif
          
       endif else begin

          mould[iobj].v_infall = model.v_fid*(cos(theta) + $
            (mould[iobj].v_LG - model.v_LG*cos(theta))/r0a*(r0a/model.v_LG)^(1-gamma))

          mould[iobj].v_cosmic = mould[iobj].v_LG + total(mould[iobj].v_infall)

       endelse

       z_cosmic = mould[iobj].v_cosmic/light

       mould[iobj].distance = (light/H0)*(z_cosmic-z_cosmic^2.0/2.0*(1.0+q0)) ; proper distance [Mpc]
       if not keyword_set(proper) then $
         mould[iobj].distance = mould[iobj].distance*(1+z_cosmic)             ; luminosity distance [Mpc]

    endfor

    neg = where(mould.distance lt 0.0,nneg)
    if nneg ne 0L then message, 'WARNING:  The following objects have negative distances: ', /info
    for i = 0L, nneg-1L do print, neg[i], '   ', mould[neg[i]].object, mould[neg[i]].distance

return, mould
end
