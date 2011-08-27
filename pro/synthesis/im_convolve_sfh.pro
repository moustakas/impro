;+
; NAME:
;   IM_CONVOLVE_SFH()
;
; PURPOSE:
;   Convolve an SSP with an arbitrary star formation history.
;
; INPUTS: 
;   ssp - structure describing the SSP grid:
;     age  - age vector (yr) [NSSP]
;     wave - wavelength vector [NPIX]
;     flux - flux vector [NPIX,NSSP]
;
; OPTIONAL INPUTS:
;   tau - characteristic time for exponentially declining star
;     formation (Gyr) (stronger than SFH and TIME)
;   mstar - stellar mass of the SSP with time [NSSP]
;   nsamp - time oversampling factor (default 2)
; 
; OPTIONAL INPUTS/OUTPUTS (see COMMENTS):
;   sfh - star formation rate at each TIME (Msun/yr) [NSFH]
;   time - time vector corresponding to SFH (Gyr) [NSFH] 
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
;   The star formation history can be specified using *either* TAU,
;   *or* SFH and TIME.  If TAU is specified then we assume a standard
;   star formation history of the form SFR=exp(-TIME/TAU)/TAU
;   (corresponding to one solar mass of stars formed after infinite
;   time).  
;
;   Alternatively, an arbitrary star formation history can be passed
;   using SFH and TIME.  Note that the input history should have good
;   time resolution during rapid variations in the star formation
;   rate, otherwise significant inaccuracies could occur.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 20, UCSD - written, with input from
;     C. Tremonti and A. Diamond-Stanic
;   jm11mar04ucsd - various bug fixes; better handing of the time
;     vector 
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

function im_convolve_sfh, ssp, tau=tau, sfh=sfh, time=time, $
  mstar=mstar, cspmstar=cspmstar, nsamp=nsamp, debug=debug

    if (n_elements(ssp) eq 0L) then begin
       doc_library, 'im_convolve_sfh'
       return, -1
    endif

    if (tag_indx(ssp,'age') eq -1) or (tag_indx(ssp,'wave') eq -1) or $
      (tag_indx(ssp,'flux') eq -1) then begin
       splog, 'Incompatible SSP structure format'
       return, -1
    endif
    npix = n_elements(ssp.wave)

; build the SFH
    nsfh = n_elements(sfh)
    ntau = n_elements(tau)
    if (nsfh eq 0) and (ntau eq 0) then begin
       splog, 'Must specify either TAU *or* SFH and TIME!'
       return, -1
    endif
    if (ntau ne 0) then begin ; stronger than SFH and TIME
       if (ntau gt 1) then begin
          splog, 'TAU must be a scalar'
          return, -1
       endif
       time = ssp.age/1D9
       if (tau eq 0.0) then sfh = time*0.0 else $
         sfh = exp(-time/tau)/(tau*1D9) ; [Msun/yr]
       nsfh = n_elements(sfh)
    endif else begin
       if (nsfh eq 0) or (n_elements(time) eq 0) then begin
          splog, 'SFH and TIME inputs required'
          return, -1
       endif
       if (nsfh ne n_elements(time)) then begin
          splog, 'Incompatible dimensions of SFH and TIME'
          return, -1
       endif
    endelse

    if (n_elements(nsamp) eq 0) then nsamp = 2 else nsamp = fix(nsamp>2)

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

; special case of TAU=0
    if (ntau gt 0) then begin
       if (tau eq 0.0) then begin
          if (nmstar ne 0) then cspmstar = mstar
          cspflux = ssp.flux
          return, cspflux
       endif
    endif
    
; integrate over age
    if (ntau eq 0) then begin
       bigtime = [ssp.age,time*1D9]
       bigtime = bigtime[uniq(bigtime,sort(bigtime))]
    endif else bigtime = ssp.age

;   t0 = systime(1)
    for ii = 0, nsfh-1 do begin
; build the oversampled time array
       age = time[ii]*1D9 ; [yr]
       otime = [0D,bigtime<age,(age-bigtime)>0,age]
       otime = otime[uniq(otime,sort(otime))]
       thistime = interpolate(otime,dindgen(nsamp*n_elements(otime)-(nsamp-1))/(nsamp*1D))
       nthistime = n_elements(thistime)
       
; interpolate       
       sspindx = findex(ssp.age,reverse(thistime))
       isspflux = interpolate(ssp.flux,sspindx,/grid)

       if (ntau gt 0) then thissfh = exp(-thistime/(tau*1D9))/(tau*1D9) else $ ; [Msun/yr]
         thissfh = interpolate(sfh,findex(time*1D9,thistime))                  ; [Msun/yr]
       if keyword_set(debug) then begin
          if (ii eq 0) then djs_plot, time, sfh, psym=6, xsty=3, ysty=3, /xlog
          djs_oplot, thistime, thissfh, psym=6, sym=0.2, color='orange'
;         cc = get_kbrd(1)
       endif

; do the convolution integral       
       dt = (shift(thistime,-1)-shift(thistime,+1))/2.0
       dt[0] = (thistime[1]-thistime[0])/2.0
       dt[nthistime-1] = (thistime[nthistime-1]-thistime[nthistime-2])/2.0

       weight = thissfh*dt ; [Msun]
       vweight = rebin(reform(weight,1,nthistime),npix,nthistime) 
       cspflux[*,ii] = total(isspflux*vweight,2,/double)

; compute the output stellar mass
       if (nmstar gt 0) then cspmstar[ii] = $
         total(interpolate(mstar,sspindx)*weight,/double)
    endfor
;   splog, 'Time = ', systime(1)-t0

return, cspflux
end
