;+
; NAME:
;   LF_FIT_DOUBLE_POWERLAW()
;
; PURPOSE:
;   Fit a binned luminosity function using a double power-law model. 
;
; INPUTS: 
;   loglum - log-base-10 luminosity at the center of each bin [NGAL] 
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
;   powerlaw - output structure with
;     .PHISTAR - number density at the inflection point of the LF 
;     .LOGLSTAR - luminosity at the inflection point of the LF 
;     .ALPHA - first power law slope
;     .BETA - second power law slope
;     .PHISTAR_ERR
;     .LOGLSTAR_ERR
;     .ALPHA_ERR
;     .BETA_ERR
;     .CHI2_DOF - reduced chi^2
;     .COVAR[4,4] - covariance matrix
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Apr 14, NYU
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

function lf_fit_powerlaw_double_func, x, params
return, lf_double_powerlaw(x,params[0],params[1],params[2],params[3])
end

function lf_fit_double_powerlaw, loglum, phi, phierr, parinfo=parinfo, quiet=quiet

    ngal = n_elements(loglum)
    if (ngal eq 0L) or (ngal ne n_elements(phi)) or $
      (ngal ne n_elements(phierr)) then begin
       doc_library, 'lf_fit_powerlaw'
       return, -1
    endif
    
; initialize the parameter structure
    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, fixed: 0, limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,4)
       parinfo[0].value = 1D-2
       parinfo[1].value = alog10(5D9)
       parinfo[2].value = +0.5D
       parinfo[3].value = +2D
    endif

; do the fit    
    params = mpfitfun('lf_fit_powerlaw_double_func',loglum,phi,$
      phierr,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof,covar=covar)

; pack it in and return    
    powerlaw = {phistar: params[0], loglstar: params[1], $
      alpha: params[2], beta: params[3], phistar_err: perror[0], $
      loglstar_err: perror[1], alpha_err: perror[2], $
      beta_err: perror[3], covar: covar, chi2_dof: chi2/(dof+(dof eq 0))}

return, powerlaw
end
