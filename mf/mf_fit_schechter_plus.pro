;+
; NAME:
;   mf_fit_schechter_plus
; PURPOSE:
;   fit "double" schechter function to set of points
; USAGE:
;   mf_fit_schechter_plus, mass, phi, phierr, schechter_plus [, mden=]
; INPUTS:
;   mass - stellar mass at center of each bin
;   phi - luminosity function at each bin
;   phierr - error in l.f. at each bin
; OUTPUTS:
;   schechter - structure with
;                  .MSTAR 
;                  .PHISTAR 
;                  .ALPHA 
;                  .PHIPLUS 
;                  .ALPHAPLUS 
;                  .MSTAR_ERR
;                  .PHISTAR_ERR
;                  .ALPHA_ERR
;                  .PHIPLUS_ERR
;                  .ALPHAPLUS_ERR
; OPTIONAL OUTPUTS:
;   mden - total mass density in based on the best fit
; COMMENTS:
;   The version of the double schechter function we use here is the
;   sum of two Schecter functions with the same exponential cutoff but
;   different normalizations and faint-end slopes. Returns the steeper
;   faint end slope term in PHIPLUS and ALPHAPLUS. 
; REVISION HISTORY:
;   2010-Feb-11 J. Moustakas, UCSD, largely based on
;     M. Blanton's LF_FIT_SCHECHTER_PLUS 
;-

function mf_fit_schechter_plus_func, x, params
return, mf_schechter_plus(x,params[0],params[1],params[2],$
  params[3],params[4])
end

pro mf_fit_schechter_plus, logmass, phi, phistddev, $
  schechter_plus, parinfo=parinfo, quiet=quiet

    if(n_tags(schechter_plus) eq 0) then begin
       schechter_plus = {phistar:1.0D-2, mstar: 10.5D, alpha:-1.0D, phiplus:1.0D-3, $
         alphaplus:-1.0D, phistar_err: 0.0D, mstar_err: 0.0D, $
         alpha_err: 0.0D, phiplus_err: 0.0D, alphaplus_err: 0.0D, $
         chi2_dof: 1E6, rho: 0.0D, rho_err: 0.0D}
    endif

    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,5)
       parinfo[0].value = schechter.phistar
       parinfo[0].limited[0] = 1
       parinfo[0].limits[0] = 1D-8

       parinfo[3].value = schechter.phiplus
       parinfo[3].limited[0] = 1
       parinfo[3].limits[0] = 1D-8

       parinfo[1].value = schechter.mstar
       parinfo[2].value = schechter.alpha
       parinfo[4].value = schechter.alphaplus
    endif

    params = mpfitfun('mf_fit_schechter_plus_func',logmass,phi,$
      phistddev,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof)
    if (dof gt 0.0) then schechter_plus.chi2_dof = chi2/float(dof)

    schechter_plus.phistar = params[0]
    schechter_plus.mstar = params[1]
    schechter_plus.alpha = params[2]
    schechter_plus.phiplus = params[3]
    schechter_plus.alphaplus = params[4]
    schechter_plus.phistar_err = perror[0]
    schechter_plus.mstar_err = perror[1]
    schechter_plus.alpha_err = perror[2]
    schechter_plus.phiplus_err = perror[3]
    schechter_plus.alphaplus_err = perror[4]

return
end
