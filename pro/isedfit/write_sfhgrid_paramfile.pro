;+
; NAME:
;   WRITE_SFHGRID_PARAMFILE()
;
; PURPOSE:
;   Initialize the global parameter file for iSEDfit.
;
; INPUTS: 
;   sfhgrid_paramfile - output parameter file name (e.g.,
;     PREFIX+'_sfhgrid.par')
;
; OPTIONAL INPUTS: 
;   zlog - distribute redshifts logarithmically in the range
;     [MINZ,MAXZ]? [0=no, 1=yes] (default 0)
;   h100 - Hubble constant relative to 100 km/s/Mpc (default 0.7)
;   omega0 - matter density (default 0.3)
;   omegal - vacuum energy density (default 0.7)
;   igm - include IGM attenuation in all calculations? [0=no, 1=yes]
;     (default 1)
;   isedpath - full path name to the input and output files (default ./) 
;
; KEYWORD PARAMETERS:
;   clobber - overwrite existing parameter file
; 
; OUTPUTS: 
;   This code writes a parameter file called SFHGRID_PARAMETER.PAR and
;   will also optionally return the PARAMS structure. 
; 
; COMMENTS:
;   IGM attenuation should only be necessary if the bluest filters
;   under consideration lie blueward of 1215 A at the highest
;   redshifts.  Note that some calculations are slower if IGM
;   attenuation is included.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2012 Sep 19, Siena 
;
; Copyright (C) 2012, John Moustakas
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

function init_sfhgrid
; assign default no-burst parameters
    params = {$
      sfhgrid:          -1L,$
      nage:             50L,$
      nmonte:         1000L,$
      tau:       [0.01,1.0],$
      Z:       [0.004,0.04],$
      AV:        [0.35,2.0],$ ; gamma-distribution parameters
      mu:           [0.1,4],$ ; gamma-distribution parameters
      fburst:    [0.03,4.0],$ ; 
      dtburst:   [0.03,0.3],$ ; [Gyr]
      trunctau: [-1.0,-1.0],$
      minage:           0.1,$ ; [Gyr]
      maxage:          13.0,$ ; [Gyr]
      mintburst:       -1.0,$
      maxtburst:       -1.0,$
      pburst:          -1.0,$
      interval_pburst:  2.0,$ ; [Gyr]
      fractrunc:       -1.0,$
      oneovertau:        1L,$ ; note!
      delayed:          -1L,$
      flatAV:           -1L,$
      flatmu:           -1L,$
      flatfburst:       -1L,$
      flatdtburst:      -1L,$
      bursttype:         1L}  ; Gaussian burst
return, params
end    

pro write_sfhgrid_paramfile, sfhgrid_paramfile, params, sfhgrid=sfhgrid, nage=nage, $
  nmonte=nmonte, tau=tau, Z=Z, AV=AV, mu=mu, fburst=fburst, dtburst=dtburst, $
  trunctau=trunctau, minage=minage, maxage=maxage, mintburst=mintburst, $
  maxtburst=maxtburst, pburst=pburst, interval_pburst=interval_pburst, $
  fractrunc=fractrunc, oneovertau=oneovertau, delayed=delayed, flatAV=flatAV, $
  flatmu=flatmu, flatfburst=flatfburst, flatdtburst=flatdtburst, bursttype=bursttype, $
  append=append, clobber=clobber, preset_bursts=preset_bursts

    if n_elements(sfhgrid_paramfile) eq 0 then begin
       doc_library, 'write_sfhgrid_paramfile'
       return
    endif
    
; build the parameter structure, which by default does not have bursts 
    params = init_sfhgrid()

; optionally add bursts
    if keyword_set(preset_bursts) then params.pburst = 0.5

; override any of the free parameters    
    if n_elements(nage) ne 0 then params.nage = nage
    if n_elements(nmonte) ne 0 then params.nmonte = nmonte
    if n_elements(mintburst) ne 0 then params.mintburst = mintburst
    if n_elements(maxtburst) ne 0 then params.maxtburst = maxtburst
    if n_elements(pburst) ne 0 then params.pburst = pburst
    if n_elements(interval_pburst) ne 0 then params.interval_pburst = interval_pburst
    if n_elements(fractrunc) ne 0 then params.fractrunc = fractrunc
    if n_elements(bursttype) ne 0 then params.bursttype = bursttype

    if keyword_set(delayed) then begin
       params.delayed = 1
       params.oneovertau = 0
    endif
    params.oneovertau = keyword_set(oneovertau)
    params.flatAV = keyword_set(flatAV)
    params.flatmu = keyword_set(flatmu)
    params.flatfburst = keyword_set(flatfburst)
    params.flatdtburst = keyword_set(flatdtburst)

    if n_elements(minage) ne 0 then begin
       params.minage = minage
       if n_elements(mintburst) eq 0 and params.pburst gt 0 then params.mintburst = minage
    endif
    if n_elements(maxage) ne 0 then begin
       params.maxage = maxage
       if n_elements(maxtburst) eq 0 and params.pburst gt 0 then params.maxtburst = maxage
    endif

    if params.delayed and params.oneovertau then begin
       splog, 'DELAYED and ONEOVERTAU may not work well together; choose one!'
       return
    endif

    if n_elements(tau) ne 0 then begin
       if n_elements(tau) ne 2 then message, 'TAU must be a 2-element array!'
       params.tau = tau
    endif
    if n_elements(AV) ne 0 then begin
       if n_elements(AV) ne 2 then message, 'AV must be a 2-element array!'
       params.AV = AV
    endif
    if n_elements(Z) ne 0 then begin
       if n_elements(Z) ne 2 then message, 'Z must be a 2-element array!'
       params.Z = Z
    endif
    if n_elements(mu) ne 0 then begin
       if n_elements(mu) ne 2 then message, 'MU must be a 2-element array!'
       params.mu = mu
    endif
    if n_elements(fburst) ne 0 then begin
       if n_elements(fburst) ne 2 then message, 'FBURST must be a 2-element array!'
       params.fburst = fburst
    endif
    if n_elements(dtburst) ne 0 then begin
       if n_elements(dtburst) ne 2 then message, 'DTBURST must be a 2-element array!'
       params.dtburst = dtburst
    endif
    if n_elements(trunctau) ne 0 then begin
       if n_elements(trunctau) ne 2 then message, 'TRUNCTAU must be a 2-element array!'
       params.trunctau = trunctau
    endif
    
; overwrite or append?  assign unique SFHGRID numbers
    if keyword_set(append) then begin
       if file_test(sfhgrid_paramfile) eq 0 then begin
          splog, 'Parameter file '+sfhgrid_paramfile+' does not exist; unable to APPEND.'
          return
       endif
       splog, 'Appending to '+sfhgrid_paramfile
       params1 = yanny_readone(sfhgrid_paramfile)
       if n_elements(sfhgrid) eq 0 then params.sfhgrid = max(params1.sfhgrid)+1 else $
         params.sfhgrid = sfhgrid
       params = [params1,params]
    endif else begin
       if n_elements(sfhgrid) eq 0 then params.sfhgrid = 1 else $
         params.sfhgrid = sfhgrid
    endelse

    uu = uniq(params.sfhgrid,sort(params.sfhgrid))
    if n_elements(uu) ne n_elements(params) then message, $
      'SFHGRID numbers must be unique'
    
    if im_file_test(sfhgrid_paramfile,clobber=keyword_set(clobber) or $
      keyword_set(append)) then return

; write out
    hdr = [['# iSEDfit star formation history grid (SFHGRID) parameters'],$
      ['# Generated by WRITE_SFHGRID_PARAMFILE on '+im_today()]]
    splog, 'Writing '+sfhgrid_paramfile
    yanny_write, sfhgrid_paramfile, ptr_new(params), $
      stnames='SFHGRIDPARAMS', /align, hdr=hdr

return
end
