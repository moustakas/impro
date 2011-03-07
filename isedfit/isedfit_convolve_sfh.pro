;+
; NAME:
;   ISEDFIT_CONVOLVE_SFH()
;
; PURPOSE:
;   Convolve an SSP with an arbitrary star formation history. 
;
; INPUTS: 
;   ssp - structure describing the SSP grid:
;     age  - age vector in years [NAGE]
;     wave - wavelength vector [NPIX]
;     flux - flux vector [NPIX,NAGE]
;
; INPUTS:
;   infosfh - iSEDfit-style star formation history structure
;     tau - characteristic e-folding time [Gyr]
;     nburst - number of bursts
;     tburst - time each burst begins [Gyr]
;     dtburst - burst duration [Gyr]
;     fburst - mass fraction of each burst
;     tauburst - truncated burst timescale [Gyr]
;     ebv, mu - (optional) reddening factors
;
; OPTIONAL INPUTS:
;   time - desired output age vector; alternatively, can pass NSFH,
;     MINTIME, and MAXTIME, and TIME will contain the age vector
;     corresponding to SFH
;   nsfh, mintime, maxtime - number of elements in the output SFH and
;     TIME vectors (default 100), as well as the minimum and maximum
;     time/age (default to the minimum and maximum ages of the input
;     SSP)
;   mstar - stellar mass of the SSP with time [NAGE]
;   tbc - birth-cloud dissipation timescale, if using the Charlot &
;     Fall (2000) attenuation curve (default 10 Myr)
;
; OPTIONAL OUTPUTS:
;   sfh - star formation rate at each TIME [NSFH]
;   cspmstar - stellar mass of the CSP with time [NSFH]
;   calzetti, odonnell, charlot, smc - attenuate the output spectrum
;     by the amount of reddening in the 
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
;   J. Moustakas, 2011 Jan 27, UCSD
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

function isedfit_convolve_sfh, ssp, infosfh=infosfh, time=time, $
  sfh=sfh, nsfh=nsfh, mintime=mintime, maxtime=maxtime, mstar=mstar, $
  cspmstar=cspmstar, nsamp=nsamp, tbc=tbc, charlot=charlot, $
  odonnell=odonnell, calzetti=calzetti, smc=smc

    if (n_elements(ssp) eq 0L) then begin
       doc_library, 'isedfit_convolve_sfh'
       return, -1
    endif

    if (tag_indx(ssp,'age') eq -1) or (tag_indx(ssp,'wave') eq -1) or $
      (tag_indx(ssp,'flux') eq -1) then begin
       splog, 'Incompatible SSP structure format'
       return, -1
    endif

; check for the SFH    
    if (n_elements(time) eq 0) then begin
       if (n_elements(nsfh) eq 0) then nsfh = 100
       if (n_elements(mintime) eq 0) then mintime = min(ssp.age)/1D9
       if (n_elements(maxtime) eq 0) then maxtime = max(ssp.age)/1D9
       time = build_isedfit_agegrid(infosfh,minage=mintime,$
         maxage=maxtime,nage=nsfh)
    endif else nsfh = n_elements(time)
    if (n_elements(sfh) eq 0) then sfh = isedfit_reconstruct_sfh(infosfh,age=time)
    if (n_elements(sfh) ne n_elements(time)) then begin
       splog, 'Incompatible dimensions of SFH and TIME!'
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
    if (n_elements(nsamp) eq 0) then nsamp = 2.0
    nthistime = long((nsamp*(1+findgen(nsfh)))>10)

; deal with reddening, if requested
    if (n_elements(tbc) eq 0) then tbc = 0.01D ; [10 Myr]
    if (keyword_set(charlot) eq 0) and (keyword_set(odonnell) eq 0) and $
      (keyword_set(calzetti) eq 0) and (keyword_set(smc) eq 0) then calzetti = 1 ; default
    kl = k_lambda(ssp.wave,charlot=charlot,odonnell=odonnell,$
      calzetti=calzetti,smc=smc,/silent) ; attenuation curve
    kl = reform(kl,npix,1)

; deal with the first time bin, including possible time-dependent
; attenuation
    sspindx = findex(ssp.age,time[0]*1D9)
    weight = 1D9*time[0]*sfh[0] ; [Msun]

    if keyword_set(charlot) then begin
       if (time[0] le tbc) then dust = 10^(-0.4*infosfh.ebv*kl) else $
         dust = 10^(-0.4*infosfh.mu*infosfh.ebv*kl)
    endif else dust = 1.0
    
    cspflux[*,0] = replicate(weight,npix)*dust*interpolate(ssp.flux,lindgen(npix),sspindx,/grid)
    if (nmstar gt 0) then cspmstar[0] = interpolate(mstar,sspindx)*weight
    if (infosfh.tau eq 0) then begin ; special case
       cspflux[*,0] += dust*interpolate(ssp.flux,lindgen(npix),sspindx,/grid)
       if (nmstar gt 0) then cspmstar[0] += interpolate(mstar,sspindx)
    endif

; unfortunately, we have to loop over time!
;   t0 = systime(1)
    for ii = 1L, nsfh-1 do begin
       age = [0D,build_isedfit_agegrid(infosfh,minage=min(time),$
         maxage=time[ii],nage=nthistime[ii]-1,/linear)] ; linear!
       thistime = abs(reverse(time[ii]-(age-min(age)))) ; avoid negative roundoff
if ii eq 132 then stop       
       iindx = findex(time,thistime)
       thissfh = interpolate(sfh,iindx) ; [Msun/yr]

       sspindx = findex(ssp.age,reverse(age)*1D9)
       isspflux = interpolate(ssp.flux,lindgen(npix),sspindx,/grid)

;      djs_plot, thistime, thissfh, psym=-6, /xlog, xrange=[min(time),time[ii]]
;      cc = get_kbrd(1)
       
; do the convolution integral       
       dt = (shift(thistime,-1)-shift(thistime,+1))/2.0
       dt[0] = (thistime[1]-thistime[0])/2.0
       dt[nthistime[ii]-1] = (thistime[nthistime[ii]-1]-$
         thistime[nthistime[ii]-2])/2.0

       weight = thissfh*dt*1D9 ; [Msun]
       vweight = rebin(reform(weight,1,nthistime[ii]),npix,nthistime[ii])

; time-dependent attenuation?
       if keyword_set(charlot) then begin
          young = where(age le tbc,nyoung,comp=old,ncomp=nold)
          ebv = fltarr(npix,nthistime[ii])
          if (nyoung ne 0) then ebv[*,young] = infosfh.ebv
          if (nold ne 0) then ebv[*,old] = infosfh.mu*infosfh.ebv
          dust = 10^(-0.4*ebv*rebin(kl,npix,nthistime[ii]))
;         if ii eq 30 then stop
       endif else dust = 1.0

; add it up!       
       cspflux[*,ii] = total(isspflux*vweight*dust,2,/double)
       if (nmstar gt 0) then cspmstar[ii] = total(interpolate(mstar,sspindx)*weight,/double)

       if (infosfh.tau eq 0) then begin ; special case
          if keyword_set(charlot) then begin
             if (max(age) le tbc) then dust = 10^(-0.4*infosfh.ebv*kl) else $
               dust = 10^(-0.4*infosfh.mu*infosfh.ebv*kl)
          endif else dust = 1.0
          cspflux[*,ii] += dust*interpolate(ssp.flux,lindgen(npix),sspindx[0],/grid)
          if (nmstar gt 0) then cspmstar[ii] += interpolate(mstar,sspindx[0])
       endif
;      if ii eq 60 then stop
    endfor
;   splog, systime(1)-t0

; final attenuation
    if (keyword_set(charlot) eq 0) then begin
       dust = 10^(-0.4*infosfh.ebv*rebin(kl,npix,nsfh))
       cspflux = cspflux*dust
    endif
    
return, cspflux
end
