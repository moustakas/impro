;+
; NAME:
;   MF_FIT_MODIFIED_SCHECHTER()
;
; PURPOSE:
;   Fit a binned stellar mass function using a modified Schechter
;   model.
;
; INPUTS: 
;   logmass - log-base-10 stellar mass at the center of each bin
;     (presumably in units of Msun) [NGAL] 
;   phi - number density, i.e., the mass function (Mpc^-3) [NGAL] 
;   phierr - 1-sigma uncertainty in PHI in the same units [NGAL]
;
; OPTIONAL INPUTS: 
;   parinfo - MPFIT() parameter structure; can be used to constrain or
;     fix any of the model parameters.
;
; KEYWORD PARAMETERS: 
;   quiet - suppress MPFIT() messages
;
; OUTPUTS: 
;   modified_schechter - output structure with
;     .PHISTAR - number density at the 'knee' of the MF
;     .LOGMSTAR - log-base-10 stellar mass at the 'knee' of the MF
;     .ALPHA - low-mass slope
;     .BETA - second power-law parameter
;     .PHISTAR_ERR
;     .LOGMSTAR_ERR
;     .ALPHA_ERR
;     .BETA_ERR
;     .CHI2_DOF - reduce chi^2
;     .COVAR[3,3] - covariance matrix
;
; COMMENTS:
;   The model is the sum of two Schecter functions with the same
;   exponential cutoff but different normalizations and faint-end slopes.
;
;   MPFIT() has difficulty getting the errors right if MSTAR
;   isn't logarithmic. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 11, UCSD - based largely on
;     M. Blanton's LF_FIT_MODIFIED_SCHECHTER
;
; Copyright (C) 2010, John Moustakas
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

function mf_fit_modified_schechter_func, x, params
return, mf_modified_schechter(x,params[0],params[1],params[2],params[3])
end

function mf_fit_modified_schechter, logmass, phi, phierr, parinfo=parinfo, quiet=quiet

    ngal = n_elements(logmass)
    if (ngal eq 0L) or (ngal ne n_elements(phi)) or $
      (ngal ne n_elements(phierr)) then begin
       doc_library, 'mf_fit_modified_schechter'
       return, -1
    endif
    
; initialize the parameter structure
    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, limited: [0,0], limits: [0.0D,0.0D]}

       parinfo = replicate(parinfo,5)
       parinfo[0].value = 1D-2 ; phi*
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-5

       parinfo[3].value = 1D-3 ; phi+
       parinfo[3].limited[0] = 1
       parinfo[3].limits[0] = 1D-5

       parinfo[1].value = 10.5D ; log(M*)
       parinfo[1].limited = 1
       parinfo[1].limits = [8D,12.5D]

       parinfo[2].value = -0.9D ; alpha
       parinfo[2].limited[1] = 1
       parinfo[2].limits[0] = 0D

       parinfo[4].value = -1.5D ; alpha+
       parinfo[4].limited[1] = 1
       parinfo[4].limits[0] = 0D
    endif

; do the fit    
    params = mpfitfun('mf_fit_modified_schechter_func',logmass,phi,$
      phierr,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof,covar=covar)

    modified_schechter = {phistar: params[0], logmstar: params[1], $
      alpha: params[2], beta: params[3], phistar_err: perror[0], $
      logmstar_err: perror[1], alpha_err: perror[2], $
      beta_err: perror[3], covar: covar, chi2_dof: chi2/(dof+(dof eq 0))}

return, modified_schechter
end
