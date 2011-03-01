;+
; NAME:
;       IM_HIIREGION_DEPROJECT()
;
; PURPOSE:
;       Compute the de-projected galactocentric position of an HII
;       region given its nuclear-relative offset coordinates and the
;       inclination and position angles of the parent galaxy.  The
;       position angle is measured positive from North to East.
;
; INPUTS:
;       incl  - galaxy inclination angle [degree]
;       pa    - galaxy position angle [N-->E, degree]
;       raoff - HII region right ascension offset [arcsec]
;       deoff - HII region declination offset [arcsec]
;
; OPTIONAL INPUTS:
;       hii_region - optional HII region name
;
; KEYWORD PARAMETERS:
;       help - print the syntax of this routine
;
; OUTPUTS:
;       hii_radius - de-projected galactocentric radius [arcsec] 
;
; OPTIONAL OUTPUTS:
;       hii_phi - de-projected HII-region position angle [degree] 
;
; COMMENTS:
;       This isn't very smart, but if RAOFF *and* DEOFF are equal to
;       -999, then it is assumed that the coordinates are not known.
;       In principal, these quantities can be less than -999 (e.g.,
;       the HII regions in M31 from Bresolin, Kennicutt, & Garnett
;       1999).  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Nov 08, U of A, originally written based on
;          a FORTRAN routine by F. Bresolin. 
;       jm06jan19uofa - documented, additional error checking,
;          streamlined 
;       jm06jan23uofa - totally re-written and documentation updated  
;
; Copyright (C) 2005-2006, John Moustakas
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

function im_hiiregion_deproject, incl1, pa1, raoff, deoff, $
  hii_phi=hii_phi, hii_region=hii_region, help=help

    if keyword_set(help) then begin
       print, 'Syntax - hii_radius = im_hiiregion_deproject(incl,pa,$'
       print, '   raoff,deoff,hii_phi=,hii_region=,/help)'
       return, -1L
    endif

    nobject = n_elements(raoff)
    if (n_elements(deoff) ne nobject) then begin
       print, 'RAOFF and DEOFF have incompatible dimensions!'
       return, -1L
    endif

    case n_elements(incl1) of
       0L: begin
          print, 'INCL input required:'
          junk = im_hiiregion_deproject(/help)
          return, -1L
       end
       1L: incl = replicate(incl1,nobject)
       else: begin
          if (n_elements(incl1) ne nobject) then begin
             print, 'RAOFF, DEOFF, and INCL have incompatible dimensions!'
             return, -1L
          endif else incl = incl1
       end
    endcase
    
    case n_elements(pa1) of
       0L: begin
          print, 'PA input required:'
          junk = im_hiiregion_deproject(/help)
          return, -1L
       end
       1L: pa = replicate(pa1,nobject)
       else: begin
          if (n_elements(pa1) ne nobject) then begin
             print, 'RAOFF, DEOFF, and PA have incompatible dimensions!'
             return, -1L
          endif else pa = pa1
       end
    endcase
    
    if (n_elements(hii_region) eq 0L) then hii_region = 'Region'+$
      string(lindgen(nobject),format='(I0)')

; call this routine recursively    
    
    if (nobject gt 1L) then begin
       hii_radius = fltarr(nobject)
       hii_phi = hii_radius*0.0
       for i = 0L, nobject-1L do begin
          hii_radius[i] = im_hiiregion_deproject(incl[i],pa[i],raoff[i],$
            deoff[i],hii_phi=hii_phi1,hii_region=hii_region[i])
          hii_phi[i] = hii_phi1
;         if hii_radius[i] gt -900.0 then print, i, $
;           hii_region[i], hii_radius[i], hii_phi[i]
       endfor
       return, hii_radius
    endif

; check for insufficient information    
    
    if (incl lt -900.0) or (pa lt -900.0) or $
      ((raoff eq -999.0) and (deoff eq -999.0)) then begin
       hii_phi = -999.0
;      print, 'insufficient information to de-project!'
       return, -999.0
    endif
    
; check for inclination>89
    
    if (incl ge 89.0) then begin
       print, 'warning: i>89 for object '+hii_region+'.'
       hii_phi = -999.0
       return, -999.0
    endif
    
; rotate the coordinate system such that the y' axis is along the
; major axis of the galaxy; deproject the galaxy disk onto the new
; y'', x''=x'; and finally compute the galactocentric distance in
; arcseconds

    x1 = -raoff ; east is usually positive; switch the coordinate axis
    y1 =  deoff

; let theta be the angle relative to the x-axis between the cardinal
; (unrotated) xy coordinate axes, and x'y' be the rotated coordinate
; axes (e.g., arfken & webber, pg. 8); we have
;
;   x' =  x*cos(theta) + y*sin(theta)
;   y' = -x*sin(theta) + y*cos(theta)
;
; since astronomical position angles are measured relative to the
; y-axis, we also have theta = pa - 90

    theta = (pa[0] - 90.0)*!dtor
    xp =  x1*cos(theta) + y1*sin(theta)
    yp = -x1*sin(theta) + y1*cos(theta)

    ypp = yp/cos(incl[0]*!dtor)   ; de-project
    hii_radius = sqrt(xp^2+ypp^2) ; [arcsec]

    if (x1 eq 0.0) then hii_phi = 0.0 else hii_phi = atan(ypp/xp)*!radeg

return, hii_radius
end
    
