;+
; NAME:
;       FLOW_DISTANCE()
;
; PURPOSE:
;       Return Virgo-centric infall galaxy distances based on both the
;       Tully (1988) calculations and the Kraan-Korteweg (1986)
;       calculations. 
;
; CALLING SEQUENCE:
;       result = flow_distance(ra,dec,H0=,searchradius=,$
;          galaxy=,/silent)
;
; INPUTS:
;       ra     - right ascension (J2000) [NOBJ]
;       dec    - declination (J2000) [NOBJ]
;
; OPTIONAL INPUTS:
;       H0           - Hubble constant (default 70) [km/s/Mpc]
;       searchradius - search radius (default: 60) [arcsec]
;       galaxy       - galaxy name (copied into RESULT)
;
; KEYWORD PARAMETERS:
;       silent - do not print messages to STDOUT
;
; OUTPUTS:
;       result - data structure with the following fields:
;          galaxy:        input galaxy name (if provided)
;          ra:            input RA
;          dec:           input DEC
;          tully_galaxy:  galaxy name in the Tully catalog
;          tully_dist:    distance in the Tully catalog [Mpc]
;          tully_match:   matching distance [arcsec]
;          kraan_galaxy:  galaxy name in the Kraan-Korteweg catalog 
;          kraan_dist:    distance in the Kraan-Korteweg catalog [Mpc] 
;          kraan_match:   matching distance [arcsec]
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       IM_HMS2DEC(), READ_88TULLY(), READ_KRAAN(), IM_DJS_ANGLE_MATCH() 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 December 10, U of A
;       jm05feb24uofa - cleaned up the code and documentation 
;
; Copyright (C) 2002-2005, John Moustakas
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

function flow_distance, ra, dec, H0=H0, searchradius=searchradius, $
  galaxy=galaxy, silent=silent

    nobj = n_elements(ra)
    if (nobj eq 0L) then begin
       print, 'Syntax - result = flow_distance(ra,dec,H0=,searchradius=,galaxy=,/silent)'
       return, -1L
    endif
    
; error checking

    if (nobj ne n_elements(dec)) then begin
       splog, 'RA and DEC have incompatible dimensions!'
       return, -1L
    endif
        
    if n_elements(searchradius) eq 0L then searchradius = 60.0 ; search radius [arcseconds]
    if not keyword_set(silent) then splog, 'Search radius (arcsec): '+$
      strtrim(string(searchradius,format='(F12.1)'),2)
    
    if (n_elements(galaxy) eq 0L) then galaxy = replicate('',nobj)

; constants
    
    if n_elements(H0) eq 0L then H0 = 70.0 ; Hubble constant [km/s/Mpc]

    virgodist = 14.6    ; HST Key Project Virgo distance (Freedman et al. 2001) [Mpc] 
    virgodist_err = 0.3 ; error [Mpc]

; initialize the output structure

    result = {$
      galaxy:         '', $
      ra:             '', $
      dec:            '', $
      tully_galaxy:   '', $
      tully_dist:  -999.0,$
      tully_match: -999.0,$
      kraan_galaxy:   '', $
      kraan_dist:  -999.0,$
      kraan_match: -999.0,$
      flow_dist:   -999.0}
    result = replicate(result,nobj)

    result.galaxy = galaxy
    result.ra = ra
    result.dec = dec

    raref = 15.0*im_hms2dec(ra)
    decref = im_hms2dec(dec)
    
; restore the Tully catalog

    tully = read_88tully()
    tra = 15.0*im_hms2dec(tully._raj2000)
    tde = im_hms2dec(tully._dej2000)

    ntot = im_djs_angle_match(raref,decref,tra,tde,dtheta=searchradius/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    good = where(mindx ne -1L,ngood)
    if not keyword_set(silent) then splog, string(ngood,format='(I3)')+'/'+$
      string(nobj,format='(I3)')+' galaxies are in Tully 1988.'

; store the Tully distances, scaling to our Hubble constant
    
    if (ngood ne 0L) then begin
       result[good].tully_galaxy = strcompress(tully[mindx[good]].name,/remove)
       result[good].tully_dist = tully[mindx[good]].r*75.0/H0 ; [Mpc]
       result[good].tully_match = mdist[good]*3600.0
    endif

; restore the Kraan-Korteweg sample of BCG's

    kraan = read_kraan()
    kraandist = virgodist*kraan.r220_1_ ; [Mpc]

    kra = 15.0*im_hms2dec(kraan._raj2000)
    kde = im_hms2dec(kraan._dej2000)
    
    ntot = im_djs_angle_match(raref,decref,kra,kde,dtheta=searchradius/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    good = where(mindx ne -1L,ngood)
    if not keyword_set(silent) then splog, 'Distances: '+string(ngood,format='(I3)')+'/'+$
      string(nobj,format='(I3)')+' galaxies are in Kraan-Korteweg.'

    if ngood ne 0L then begin
       result[good].kraan_galaxy = strcompress(kraan[mindx[good]].name1,/remove)
       result[good].kraan_dist = kraandist[mindx[good]]
       result[good].kraan_match = mdist[good]*3600.0
    endif

return, result
end    
