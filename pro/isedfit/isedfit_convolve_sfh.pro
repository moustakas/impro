;+
; NAME:
;   ISEDFIT_CONVOLVE_SFH()
;
; PURPOSE:
;   iSEDfit-specific convolution code (see IM_CONVOLVE_SFH for a
;   generalized version).
;
; INPUTS: 
;   ssp - structure describing the SSP grid:
;     age  - age vector in years [NAGE]
;     wave - wavelength vector [NPIX]
;     flux - flux vector [NPIX,NAGE]
;
;   infosfh - iSEDfit-style star formation history structure
;     tau - characteristic e-folding time [Gyr]
;     nburst - number of bursts
;     tburst - time each burst begins [Gyr]
;     dtburst - burst duration [Gyr]
;     fburst - mass fraction of each burst
;
; OPTIONAL INPUTS:
;   time - desired output age vector (Gyr)
;   mstar - stellar mass of the SSP with time [NAGE]
;   nsamp - time oversampling factor (default 2)
;
; OPTIONAL OUTPUTS:
;   sfh - star formation rate at each TIME [NSFH]
;   cspmstar - stellar mass of the CSP with time [NSFH]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS: 
;   cspflux - time-dependent spectra of the composite stellar
;     population (CSP) [NPIX,NSFH]
;
; COMMENTS:
;   Most inputs, especially TIME, should be double-precision for best
;   results. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 27, UCSD
;   jm11mar07ucsd - better age resolution; various other bug fixes
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
  sfh=sfh, mstar=mstar, cspmstar=cspmstar, nsamp=nsamp, debug=debug, $
  bigdebug=bigdebug, _extra=extra

    if (n_elements(ssp) eq 0) or (n_elements(infosfh) eq 0) then begin
       doc_library, 'isedfit_convolve_sfh'
       return, -1
    endif

    if (tag_indx(ssp,'age') eq -1) or (tag_indx(ssp,'wave') eq -1) or $
      (tag_indx(ssp,'flux') eq -1) then begin
       splog, 'Incompatible SSP structure format'
       return, -1
    endif
    npix = n_elements(ssp.wave)

    if (n_elements(nsamp) eq 0) then nsamp = 2

; check for the SFH
    if (n_elements(time) eq 0) then time = ssp.age/1D9
    sfh = isedfit_reconstruct_sfh(infosfh,outage=time,$
      debug=debug,_extra=extra)
    nsfh = n_elements(sfh)

; check for other quantities to convolve
    nmstar = n_elements(mstar)
    if (nmstar ne 0) then begin
       if (nmstar ne n_elements(ssp.age)) then begin
          splog, 'Dimensions of MSTAR and SSP.AGE do not agree'
          return, -1
       endif
       cspmstar = fltarr(nsfh)
    endif
    cspflux = fltarr(npix,nsfh)

; integrate over age
    bigtime = [ssp.age,time*1D9]
    bigtime = bigtime[uniq(bigtime,sort(bigtime))]

;   t0 = systime(1)
    for ii = 0, nsfh-1 do begin
; build the oversampled time array
       age = time[ii]*1D9 ; [yr]
       otime = [0D,bigtime<age,(age-bigtime)>0D,age]
       otime = im_double(otime) ; even though it's double, we still need this...
       otime = otime[uniq(otime,sort(otime))]
       thistime = interpolate(otime,dindgen(nsamp*n_elements(otime)-(nsamp-1))/(nsamp*1D))

       thistime = build_isedfit_agegrid(infosfh,inage=thistime/1D9)*1D9
       nthistime = n_elements(thistime)
       
; interpolate       
       sspindx = findex(ssp.age,reverse(thistime))>0
       isspflux = interpolate(ssp.flux,sspindx)
       thissfh = isedfit_reconstruct_sfh(infosfh,useage=thistime/1D9)

       if keyword_set(bigdebug) then begin
          djs_plot, time, sfh, psym=6, xsty=3, ysty=3, /xlog, /ylog, $
            yrange=[min(sfh),max(sfh)>max(thissfh)], xr=[5,10]
          djs_oplot, thistime/1D9, thissfh, psym=6, sym=0.2, color='orange'
          ww = where(max(thistime)-thistime lt 0.01D9)
          djs_oplot, thistime[ww]/1D9, thissfh[ww], psym=6, sym=0.6, color='blue'
          nb = infosfh.nburst
          if (nb gt 0) and (max(thistime/1D9) ge infosfh.tburst[0]) then begin
             t1 = findex(thistime/1D9,infosfh.tburst[0:nb-1])
             t2 = findex(thistime/1D9,infosfh.tburst[0:nb-1]+infosfh.dtburst[0:nb-1])
;            niceprint, t1, t2
;            cc = get_kbrd(1)
          endif
;         if ii eq nsfh-25 then stop          
;         print, ii, nthistime
;         cc = get_kbrd(1)
       endif

; do the convolution integral       
       dt = (shift(thistime,-1)-shift(thistime,+1))/2.0
       dt[0] = (thistime[1]-thistime[0])/2.0
       dt[nthistime-1] = (thistime[nthistime-1]-thistime[nthistime-2])/2.0

       weight = thissfh*dt ; [Msun]
       vweight = rebin(reform(weight,1,nthistime),npix,nthistime) 
       cspflux[*,ii] = total(isspflux*vweight,2)

; compute the output stellar mass
       if (nmstar gt 0) then cspmstar[ii] = $
         total(interpolate(mstar,sspindx)*weight)
    endfor
;   splog, 'Time = ', systime(1)-t0

return, cspflux
end
