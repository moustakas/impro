;+
; NAME:
;   READ_WITT_DUSTMODELS()
;
; PURPOSE:
;   Read the Witt & Gordon (2000) dust models.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   geometry - 'dusty', 'cloudy' or 'shell' (default 'shell')
;   dust_type - 'SMC' or 'MW' (default 'MW')
;   structure - 'c' (clumpy) or 'h' (homogeneous) (default 'h') 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   Structure containing all the relevant data.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2010, John Moustakas
;   J. Moustakas, 2010 Mar 17, UCSD - largely based on
;     M. Blanton's WITT_EXT()  
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

function read_witt_dustmodels, geometry=geometry, $
  dust_type=dust_type, structure=structure

    common witt_ext_common, w_lambda, w_tauv, w_tau_att, $
      w_geometry, w_dust_type, w_structure

    if (n_elements(geometry) eq 0) then geometry = 'shell'
    if (n_elements(dust_type) eq 0) then dust_type = 'MW'
    if (n_elements(structure) eq 0) then structure = 'h'
    
    if (n_elements(w_geometry) gt 0) then $
      if (w_geometry ne geometry) then w_lambda = 0
    if (n_elements(w_dust) gt 0) then $
      if (w_dust_type ne dust_type) then w_lambda = 0
    if (n_elements(w_structure) gt 0) then $
      if (w_structure ne structure) then w_lambda = 0
    if (n_elements(w_lambda) le 1) then begin
       rootdir = getenv('CATALOGS_DIR')+'/00witt/'
       k_read_ascii_table, w_lambda, rootdir+'witt.'+geometry+'_'+$
         dust_type+ '_'+structure+'.lambda.dat'
       k_read_ascii_table, w_tauv, rootdir+'witt.'+geometry+'_'+$
         dust_type+'_'+structure+'.tauv.dat'
       k_read_ascii_table, w_tau_att, rootdir+'witt.'+geometry+'_'+$
         dust_type+'_'+structure+'.tau_att.dat'
       w_geometry = geometry
       w_dust_type = dust_type
       w_structure = structure
    endif
    
; pack into a structure and return; note: A(lambda)=2.5*tau(lambda) is
; the wavelength-dependent attenuation; also, pad the attenuation
; curves on either end to ensure smooth interpolation onto, e.g., a
; galaxy SED
    dim = size(w_tau_att,/dim)
    nwave = dim[0] ; number of wavelengths
    ntauv = dim[1] ; tau_v grid

    wave = [1.0,w_lambda[0],w_lambda,w_lambda[nwave-1],1E9]
    nnwave = n_elements(wave)
    
    model = {$
      geometry:     geometry,$
      dust_type:   dust_type,$
      structure:   structure,$
      tauv:              0.0,$
      wave:             wave,$ ; [A]
      alambda: fltarr(nnwave)} ; A(lambda) [mag]
    model = replicate(model,ntauv)
    model.tauv = w_tauv
    for ii = 0, ntauv-1 do begin
       tau_att = [w_tau_att[0,ii],w_tau_att[0,ii],w_tau_att[*,ii],$
         w_tau_att[nwave-1,ii],w_tau_att[nwave-1,ii]]
       model[ii].alambda = 2.5*tau_att
       zero = where(model[ii].alambda eq 0.0,nzero) ; note!
       if (nzero ne 0) then model[ii].alambda[zero] = max(model[ii].alambda)
    endfor

return, model
end
    
