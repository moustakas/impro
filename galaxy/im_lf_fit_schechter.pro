;+
; NAME:
;   IM_LF_FIT_SCHECHTER()
;
; PURPOSE:
;   Fit a binned luminosity function using a Schechter model. 
;
; INPUTS: 
;   absmag - absolute magnitude at the center of each bin [NGAL] 
;   phi - number density, i.e., the luminosity function (Mpc^-3) [NGAL] 
;   phierr - 1-sigma uncertainty in PHI [NGAL]
;
; OPTIONAL INPUTS: 
;   parinfo - MPFIT() parameter structure; can be used to constrain or
;     fix any of the model parameters
;
; KEYWORD PARAMETERS: 
;   quiet - suppress MPFIT() messages
;
; OUTPUTS: 
;   schechter - output structure with
;     .PHISTAR - number density at the 'knee' of the LF
;     .MSTAR - absolute magnitude at the 'knee' of the LF
;     .ALPHA - low-mass slope
;     .PHISTAR_ERR
;     .MSTAR_ERR
;     .ALPHA_ERR
;     .CHI2_DOF - reduce chi^2
;     .COVAR[3,3] - covariance matrix
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Apr 14, NYU - based largely on
;     M. Blanton's LF_FIT_SCHECHTER 
;
; Copyright (C) 2009, John Moustakas
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

function im_lf_fit_schechter_func, x, params
return, im_lf_schechter(x,params[0],params[1],params[2])
end

function im_lf_fit_schechter, absmag, phi, phierr, parinfo=parinfo, quiet=quiet

    ngal = n_elements(absmag)
    if (ngal eq 0L) or (ngal ne n_elements(phi)) or $
      (ngal ne n_elements(phierr)) then begin
       doc_library, 'im_lf_fit_schechter'
       return, -1
    endif
    
; initialize the parameter structure
    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, fixed: 0, limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,3)
       parinfo[0].value = 1D-2
       parinfo[1].value = -20D
       parinfo[2].value = -1.0D
    endif

; do the fit    
    params = mpfitfun('im_lf_fit_schechter_func',absmag,phi,$
      phierr,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof,covar=covar)

; pack it in and return    
    schechter = {phistar: params[0], mstar: params[1], $
      alpha: params[2], phistar_err: perror[0], $
      mstar_err: perror[1], alpha_err: perror[2], $
      covar: covar, chi2_dof: chi2/(dof+(dof eq 0))}

return, schechter
end
