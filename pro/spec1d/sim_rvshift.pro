;+
; NAME:
;       SIM_RVSHIFT
;
; PURPOSE:
;       A simulation to test the accuracy of solving for a velocity
;       shift as a function of input spectrum S/N.
;
; CALLING SEQUENCE:
;       sim_rvshift, vtrue=, snmin=, snmax=, dsn=, minwave=, $
;          maxwave=, zans=, template=, /debug, /doplot, /silent
;
; INPUTS:
;       
; OPTIONAL INPUTS:
;       vtrue    - radial velocity shift to apply [km/s]
;       snmin    - minimum S/N (default 5)
;       snmin    - maximum S/N (default 100)
;       dsn      - S/N spacing (default 10)
;       minwave  - see IM_ZTWEAK()
;       maxwave  - see IM_ZTWEAK()
;       template - scalar template number to use
;          0L: Age = 25 Myr (BA star spectrum)
;          1L: Age = 1.4 Gyr (AF star spectrum)
;          2L: Age = 11 Gyr (GK star spectrum)
;          3L: Exp SFR, 12 Gyr, tau=5 (not supported)
;
; KEYWORD PARAMETERS:
;       debug  - debugging plot
;       doplot - generate a plot of residuals
;       silent - do not print the results to STDOUT
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       IM_ZTWEAK(), IM_READ_BC03(), GET_ELEMENT
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 14, U of A
;       jm04mar09uofa - general updates; use IM_READ_BC03() 
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

pro sim_rvshift, vtrue=vtrue, snmin=snmin, snmax=snmax, dsn=dsn, $
  minwave=minwave, maxwave=maxwave, zans=zans, template=template, $
  debug=debug, doplot=doplot, silent=silent

    light = 2.99792458D5 ; speed of light [km/s]

; generate the S/N array

    if n_elements(snmin) eq 0L then snmin = 5.0
    if n_elements(snmax) eq 0L then snmax = 60.0
    if n_elements(dsn) eq 0L then dsn = 5.0

    if not keyword_set(silent) then begin
       splog, 'Minimum S/N = '+string(snmin,format='(G0.0)')
       splog, 'Maximum S/N = '+string(snmax,format='(G0.0)')
       splog, 'S/N spacing = '+string(dsn,format='(G0.0)')
    endif
    
    sn = findgen((snmax-snmin)/dsn+1)*dsn+snmin > 5.0 ; S/N floor
    nsn = n_elements(sn)

; read the template galaxy spectrum
    
    if n_elements(template) eq 0L then template = 2L
    splog, 'Selecting template #'+string(template,format='(I0)')+'.'
    case template of
       0L: bc03 = im_read_bc03(age=0.025,minwave=3500.0,maxwave=9000.0,/silent)
       1L: bc03 = im_read_bc03(age=1.4,minwave=3500.0,maxwave=9000.0,/silent)
       2L: bc03 = im_read_bc03(age=11.0,minwave=3500.0,maxwave=9000.0,/silent)
    endcase 

    flux = bc03.flux
    wave = bc03.wave
    npix = n_elements(flux)
    
    if n_elements(vtrue) eq 0L then vtrue = +200.0 ; [km/s]

    waveshift = wave*(1+vtrue/light)
    medwave = median(waveshift)
    get_element, waveshift, medwave, windx
    
    for k = 0L, nsn-1L do begin

       print, format='("S/N #",I0,"/",I0,".",A1,$)', k+1, nsn, string(13b)

; degrade the spectrum

       ferr = flux/sqrt(flux/flux[windx])/sn[k]
       newflux = flux + randomn(seed,npix)*ferr ; perturbed flux
       ivar = 1.0/ferr^2.0

       newflux = im_gauss_broaden(wave,newflux,3.0,6.0)
       
       zans1 = im_ztweak(newflux,waveshift,flux,wave,specivar=ivar,vmin=vmin,$
         vmax=vmax,minwave=minwave,maxwave=maxwave,/silent,doplot=debug)
       if k eq 0L then zans = zans1 else zans = [ [zans], [zans1] ]
       if keyword_set(debug) then begin
          if (k eq 0L) then splog, 'Press any key to continue.'
          cc = get_kbrd(1)
       endif

    endfor 
    zans = reform(zans)

    if keyword_set(doplot) or keyword_set(debug) then begin
    
       window, 0, xs=500, ys=400
       resid = zans.vshift-vtrue
       yrange = max(abs(resid+zans.vshift_err))*[-1.2,1.2]
       ploterror, sn, resid, zans.vshift_err, ps=4, xsty=3, ysty=3, $
         xthick=2.0, ythick=2.0, charsize=1.8, charthick=2.0, xtitle='S/N', $
         ytitle='Residuals', thick=2.0, yrange=yrange
       oplot, !x.crange, [0,0], line=0, thick=2.0

    endif 
       
    if not keyword_set(silent) then struct_print, zans
    
return
end    
