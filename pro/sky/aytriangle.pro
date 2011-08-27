;+
; NAME:
;   AYTRIANGLE()
;
; PURPOSE:
;   Calculate zenith distance, hour angle, altitude, azimuth and
;   parallactic angle for an object of a given RA and DEC on a given
;   DATE at a known OBSERVATORY.    
;
; INPUTS:
;   date - date of the observations in FITS standard string
;     format, [YYYY-MM-DD] (e.g., date = '2002-02-10')
;   ra - right ascension [decimal hours]
;   dec - declination [decimal degrees]
;   obs - observatory name (e.g., 'KPNO')
;
; OUTPUTS:
;   The output data structure will contain
;     .lst - local siderial time during the night [decimal hours]
;     .ha  - hour angle [decimal hours]
;     .alt - altitude [decimal degrees]
;     .zd  - zenith distance [decimal degrees]
;     .az  - azimuth [decimal degrees]
;     .pangle - parallactic angle [decimal degrees]
;
; COMMENTS
;   See Computational Spherical Astronomy (Taff, p.13-15). 

; MODIFICATION HISTORY
;   John Moustakas, 2001 July 6, U of A
;   jm09dec09ucsd - changed to a function and cleaned up
;
; Copyright (C) 2001, 2009, John Moustakas
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

function aytriangle, date, ra, dec, obs, doplot=doplot

    if (n_params() ne 4) then begin
       doc_library, 'aytriangle'
       return, -1
    endif
    
; observatory info
    observatory, obs, obsinfo
    latitude = obsinfo.latitude

; parse the input date
    mdy = strsplit(date,'-',/extract) ; [Year,Month,Day]
    
; generate a time vector from 18:00 to 06:00 hours in 10 minute
; intervals, and then convert to LST in decimal hours
    jd_ut = timegen(6,12,start=julday(mdy[1],mdy[2],$
      mdy[0],18.0+obsinfo.tz),unit='minute',step=10)
    jd_ut = reform(jd_ut,n_elements(jd_ut))

    ct2lst, lst, (360.0-obsinfo.longitude), junk, jd_ut
    nlst = n_elements(lst)
    
; initialize the output structure
    out = {$
      lst:    0.0, $
      ha:     0.0, $
      alt:    0.0, $
      zd:     0.0, $
      az:     0.0, $
      pangle: 0.0}
    out = replicate(out,nlst)
    out.lst = lst

; hour angle (negative if the object is East of the meridian)    
    out.ha = lst - ra 
    
    dtor = !dpi/180D
    radeg = 1D/dtor
    sin_alt = sin(latitude*dtor)*sin(dec*dtor) + $
      cos(latitude*dtor)*cos(dec*dtor)*cos(15D*out.ha*dtor)
    out.alt = asin(sin_alt)*radeg ; altitude (decimal degrees)

    out.zd = 90D - out.alt ; zenith distance (decimal degrees)
    
    sin_az = cos(dec*dtor)*sin(15D*out.ha*dtor)/cos(out.alt*dtor)
    az = asin(sin_az)*radeg
    neg = where(az lt 0.0,nneg)
    if (nneg ne 0) then az[neg] = az[neg] + 360D
    out.az = az
    
    sin_pangle = cos(latitude*dtor)*sin(az*dtor)/cos(dec*dtor)
    neg = where(out.ha lt 0.0,nneg,comp=pos,ncomp=npos)
    if (nneg ne 0) then out[neg].pangle = asin(-sin_pangle[neg])
    if (npos ne 0) then out[pos].pangle = asin(sin_pangle[pos])
    out.pangle = out.pangle*radeg ; parallactic angle (decimal degrees)

return, out
end
