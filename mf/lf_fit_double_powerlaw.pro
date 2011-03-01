;+
; NAME:
;   lf_fit_double_powerlaw
; PURPOSE:
;   fit double_powerlaw function to set of points
; USAGE:
;   lf_fit_double_powerlaw, lum, phi, phierr, double_powerlaw [, mden=]
; INPUTS:
;   lum - stellar lum at center of each bin
;   phi - lum function at each bin
;   phierr - error in mf at each bin
; OUTPUTS:
;   double_powerlaw - structure with
;                  .PHISTAR 
;                  .LSTAR 
;                  .ALPHA 
;                  .BETA
;                  .PHISTAR_ERR
;                  .LSTAR_ERR
;                  .ALPHA_ERR
;                  .BETA_ERR
; OPTIONAL OUTPUTS:
;   mden - total LUM density based on the best fit
; REVISION HISTORY:
;   14-Apr-2009  Written by John Moustakas, NYU, entirely based
;     on M. Blanton's lf_fit_double_powerlaw
;-
function lf_fit_double_powerlaw_func, x, params
return, lf_double_powerlaw(x,params[0],params[1],params[2],params[3])
end
;
pro lf_fit_double_powerlaw, loglum, phi, phistddev, bestfit, $
  parinfo=parinfo, quiet=quiet
    
; remember that LSTAR is logarithmic!
    if (n_tags(bestfit) eq 0) then begin
       bestfit={phistar: 1D-2, lstar: alog10(5D9), alpha: +0.5D, $
         beta: +2.0D, phistar_err: 0D, lstar_err: 0D, alpha_err: 0D, $
         beta_err: 0D, chi2_dof: 1E6}
    endif

    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,4)
       parinfo[0].value = bestfit.phistar
       parinfo[1].value = bestfit.lstar
       parinfo[2].value = bestfit.alpha
       parinfo[3].value = bestfit.beta
    endif

    params = mpfitfun('lf_fit_double_powerlaw_func',loglum,phi,$
      phistddev,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof)
    if (dof gt 0.0) then bestfit.chi2_dof = chi2/float(dof)

    bestfit.phistar = params[0]
    bestfit.lstar = params[1]
    bestfit.alpha = params[2]
    bestfit.beta = params[3]
    bestfit.phistar_err = perror[0]
    bestfit.lstar_err = perror[1]
    bestfit.alpha_err = perror[2]
    bestfit.beta_err = perror[3]

return
end
