;+
; NAME:
;   IM_CONVOLVE_SFH()
;
; PURPOSE:
;   Convolve an SSP with an arbitrary star formation history. 
;
; INPUTS: 
;   ssp - structure describing the SSP grid:
;     age  - age vector (yr) [NAGE]
;     wave - wavelength vector [NPIX]
;     flux - flux vector [NPIX,NAGE]
;   sfh - star formation rate at each TIME (Msun/yr) [NSFH]
;   time - time vector corresponding to SFH (Gyr) [NSFH] 
;
; OPTIONAL INPUTS:
;   mstar - stellar mass of the SSP with time [NAGE]
;
; OPTIONAL OUTPUTS:
;   cspmstar - stellar mass of the CSP with time [NSFH]
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   cspflux - time-dependent spectra of the composite stellar
;     population (CSP) [NPIX,NSFH]
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 20, UCSD - written, with input from
;     C. Tremonti and A. Diamond-Stanic
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

function im_convolve_sfh, ssp, sfh=sfh, time=time, mstar=mstar, $
  cspmstar=cspmstar

    if (n_elements(ssp) eq 0L) then begin
       doc_library, 'im_convolve_sfh'
       return, -1
    endif

    if (tag_indx(ssp,'age') eq -1) or (tag_indx(ssp,'wave') eq -1) or $
      (tag_indx(ssp,'flux') eq -1) then begin
       splog, 'Incompatible SSP structure format'
       return, -1
    endif

; check for the SFH
    nsfh = n_elements(sfh)
    if (nsfh eq 0) or (n_elements(time) eq 0) then begin
       splog, 'SFH and TIME inputs required'
       return, -1
    endif
    if (nsfh ne n_elements(time)) then begin
       splog, 'Incompatible dimensions of SFH and TIME'
       return, -1
    endif

    npix = n_elements(ssp.wave)
    cspflux = fltarr(npix,nsfh)

; check for other quantities to convolve
    nmstar = n_elements(mstar)
    if (nmstar ne 0) then begin
       if (nmstar ne n_elements(ssp.age)) then begin
          splog, 'Dimensions of MSTAR and SSP.AGE do not agree'
          return, -1
       endif
       cspmstar = fltarr(nsfh)
    endif
    
; scale the number of points that we integrate with NSFH
    nsamp = 2.0
    nthistime = long((nsamp*(1+findgen(nsfh)))>10)
    
; unfortunately, we have to loop over time!
;   t0 = systime(1)
    sspindx = findex(ssp.age,time[0]*1D9)
    weight = 1D9*time[0]*sfh[0] ; [Msun]
    cspflux[*,0] = replicate(weight,npix)*interpolate(ssp.flux,sspindx,/grid)
    if (nmstar gt 0) then cspmstar[0] = interpolate(mstar,sspindx)*weight

    for ii = 1L, nsfh-1 do begin
;      age = range(min(time),time[ii],nthistime[ii],/log)
;      thistime = reverse(time[ii]-(age-min(age)))
       age = [0D,range(min(time),time[ii],nthistime[ii]-1,/log)]
       thistime = abs(reverse(time[ii]-(age-min(age)))) ; avoid negative roundoff

       iindx = findex(time,thistime)
       thissfh = interpolate(sfh,iindx) ; [Msun/yr]

       sspindx = findex(ssp.age,reverse(age)*1D9)
       isspflux = interpolate(ssp.flux,sspindx,/grid)
       
; do the convolution integral       
       dt = (shift(thistime,-1)-shift(thistime,+1))/2.0
       dt[0] = (thistime[1]-thistime[0])/2.0
       dt[nthistime[ii]-1] = (thistime[nthistime[ii]-1]-$
         thistime[nthistime[ii]-2])/2.0

       weight = thissfh*dt*1D9 ; [Msun]
       vweight = rebin(reform(weight,1,nthistime[ii]),npix,nthistime[ii])
       cspflux[*,ii] = total(isspflux*vweight,2,/double)

; compute the output stellar mass
       if (nmstar gt 0) then cspmstar[ii] = $
         total(interpolate(mstar,sspindx)*weight,/double)
    endfor
;   splog, systime(1)-t0
    
return, cspflux
end
