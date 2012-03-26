;+
; NAME:
;   MF_FIT_SCHECHTER()
;
; PURPOSE:
;   Fit a binned stellar mass function using a Schechter model. 
;
; INPUTS: 
;   logmass - log-base-10 stellar mass at the center of each bin
;     (presumably in units of Msun) [NGAL] 
;   phi - number density, i.e., the mass function (Mpc^-3) [NGAL] 
;
; OPTIONAL INPUTS: 
;   phierr - 1-sigma uncertainty in PHI in the same units [NGAL]
;   parinfo - MPFIT() parameter structure; can be used to constrain or
;     fix any of the model parameters
;
; KEYWORD PARAMETERS: 
;   quiet - suppress MPFIT() messages
;
; OUTPUTS: 
;   schechter - output structure with
;     .PHISTAR - number density at the 'knee' of the MF
;     .LOGMSTAR - log-base-10 stellar mass at the 'knee' of the MF
;     .ALPHA - low-mass slope
;     .PHISTAR_ERR
;     .LOGMSTAR_ERR
;     .ALPHA_ERR
;     .CHI2_DOF - reduce chi^2
;     .COVAR[3,3] - covariance matrix
;
; COMMENTS:
;   MPFIT() has difficulty getting the errors right if MSTAR
;   isn't logarithmic, hence the convention adopted here.  Also
;   note that PHIERR is optional, but in that case the formal MPFIT
;   errors will not be correct.
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

function mf_fit_schechter_func, x, params
return, mf_schechter(x,params[0],params[1],params[2])
end

function mf_fit_schechter, logmass, phi, phierr, parinfo=parinfo, quiet=quiet

    ngal = n_elements(logmass)
    if (ngal eq 0L) or (ngal ne n_elements(phi)) then begin
       doc_library, 'mf_fit_schechter'
       return, -1
    endif

;   if (n_elements(phierr) eq 0) then begin
;      zero = where(phi le 0.0)
;      if zero[0] ne -1 then message, 'PHI contains values that are <=0!'
;      phiweights = 1.0/phi ; uniform (Poisson) weights
;   endif else begin
;      if (n_elements(phierr) ne ngal) then message, 'PHI and PHIERR are incompatible!'
;      zero = where(phierr le 0.0)
;      if zero[0] ne -1 then message, 'PHIERR cannot be <=0!'
;      phiweights = 1.0/phierr
;   endelse
    
; initialize the parameter structure
    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, fixed: 0, limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,3)
       parinfo[0].value = 1D-2
       parinfo[1].value = 10.5D
       parinfo[2].value = -1.0D
    endif

; do the fit 
    params = mpfitfun('mf_fit_schechter_func',logmass,phi,phierr,$
      weight=phiweights,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof,covar=covar)
    if mpstatus le 0 then begin
       splog, 'Problem with the fit'
       params = parinfo.value
       perror = params*0
       dof = 1
       chi2 = 0D
       covar = dblarr(n_elements(parinfo),n_elements(parinfo))-999.0
    endif
    schechter = {phistar: params[0], logmstar: params[1], $
      alpha: params[2], phistar_err: perror[0], $
      logmstar_err: perror[1], alpha_err: perror[2], $
      covar: covar, chi2_dof: chi2/(dof+(dof eq 0)), status: mpstatus}

return, schechter
end
