;+
; NAME:
;   mf_fit_schechter
; PURPOSE:
;   fit schechter function to set of points
; USAGE:
;   mf_fit_schechter, mass, phi, phierr, schechter [, mden=]
; INPUTS:
;   mass - stellar mass at center of each bin
;   phi - mass function at each bin
;   phierr - error in mf at each bin
; OUTPUTS:
;   schechter - structure with
;                  .PHISTAR 
;                  .MSTAR 
;                  .ALPHA 
;                  .PHISTAR_ERR
;                  .MSTAR_ERR
;                  .ALPHA_ERR
; OPTIONAL OUTPUTS:
;   mden - total stellar mass density based on the best fit
; REVISION HISTORY:
;   14-Apr-2009  Written by John Moustakas, NYU, entirely based
;     on M. Blanton's lf_fit_schechter
;-

function mf_fit_schechter_func, x, params
return, mf_schechter(x,params[0],params[1],params[2])
end

pro mf_fit_schechter, logmass, phi, phierr, schechter, $
  alpha_err=alpha_err, parinfo=parinfo, quiet=quiet, $
  norhoerr=norhoerr

; remember that MSTAR is logarithmic!

; alpha_err is for a very particular application; frequently, the
; faint-end slope alpha is not well-constrained and must be fixed;
; however, when computing the uncertainty on the stellar mass density
; one would like to take the *uncertainty* on alpha into account,
; which can be accomplished using this parameter (in practice, I
; simply overwrite the output error from MPFIT, which is zero is the
; parameter is fixed)
    
    if (n_tags(schechter) eq 0) then begin
       schechter = {phistar: 1.d-2, mstar: 10.5D, alpha:-1.0d, $
         phistar_err: 0.0d, mstar_err: 0.0d, alpha_err: 0.0d, $
         chi2_dof: 1E6, rho: 0.0D, rho_err: 0.0D, $
         rho_tot: 0.0D, rho_tot_err: 0.0D}
    endif

    if (n_elements(parinfo) eq 0) then begin
       parinfo = {value: 0.0D, fixed: 0,limited: [0,0], limits: [0.0D,0.0D]}
       parinfo = replicate(parinfo,3)
       parinfo[0].value = schechter.phistar
       parinfo[1].value = schechter.mstar
       parinfo[2].value = schechter.alpha
    endif

    params = mpfitfun('mf_fit_schechter_func',logmass,phi,$
      phierr,parinfo=parinfo,perror=perror,status=mpstatus,$
      quiet=quiet,bestnorm=chi2,dof=dof,covar=covar)
    if (dof gt 0.0) then schechter.chi2_dof = chi2/float(dof)
;   splog, perror, mpstatus

    schechter.phistar = params[0]
    schechter.mstar = params[1]
    schechter.alpha = params[2]
    schechter.phistar_err = perror[0]
    schechter.mstar_err = perror[1]
    schechter.alpha_err = perror[2]

; get the total mass density by integrating the observed MF and also
; by integrating the full model from 0-->infty
    schechter.rho = im_integral(10^logmass,phi)
    schechter.rho_err = sqrt(im_integral(10^logmass,phierr^2)) ; not right!
;   schechter.rho = im_integral(10^logmass,mf_schechter(logmass,schechter))

    massaxis = range(0.0,15.0,500)
    schechter.rho_tot = im_integral(10^massaxis,$
      mf_schechter(massaxis,schechter))

; sample the full covariance matrix to get the error on rho_tot,
; accounting for the fact that alpha may have been fixed
    if parinfo[2].fixed then begin
       if (n_elements(alpha_err) eq 0) then begin
          splog, 'ALPHA_ERR is zero - errors on RHO_TOT will be underestimated!'
          use_alpha_err = 0.0
       endif else use_alpha_err = alpha_err
       use_covar = covar[0:1,0:1]
       if (schechter.phistar_err gt 0.0) and (schechter.mstar_err gt 0.0) then begin
          nrand = 1000
          rand = mrandomn(seed,use_covar,nrand)
          rand[*,0] = rand[*,0]+schechter.phistar
          rand[*,1] = rand[*,1]+schechter.mstar
          ell = covar2ellipse(use_covar,nsigma=1.0) ; get 1-sigma Schechter models

          indx = get_ellipse_indices(rand[*,0],rand[*,1],$
            major=ell.major,minor=ell.minor,angle=ell.angle,$
            xcenter=schechter.phistar,ycenter=schechter.mstar)
          if (indx[0] eq -1) then message, 'Problem here!'
; Monte Carlo
          mc_rho_tot = dblarr(n_elements(indx))
          for ij = 0L, n_elements(indx)-1 do mc_rho_tot[ij] = $
            im_integral(10^massaxis,mf_schechter(massaxis,rand[indx[ij],0],$
            rand[indx[ij],1],schechter.alpha))
          schechter.rho_tot_err = djsig(mc_rho_tot)
; debugging plot       
;         djs_plot, rand[*,0], rand[*,1], psym=6, xsty=3, ysty=3
;         djs_oplot, rand[indx,0], rand[indx,1], psym=6, color='red'
;         tvellipse, ell.major, ell.minor, schechter.phistar, /data, $
;           schechter.mstar, ell.angle, color=djs_icolor('blue')
       endif else splog, 'Errors are zero!!'
    endif else begin
;      use_alpha_err = schechter.alpha_err
       splog, 'Deal with me!'
    endelse
    
;    schechter.alpha_err = perror[2]
;    if (n_elements(alpha_err) ne 0) then schechter.alpha_err = alpha_err ; overwrite!!
;
;; get the total stellar mass density; get the uncertainty using a
;; simple Monte Carlo method
;    if (keyword_set(norhoerr) eq 0) then begin
;       nmonte = 500
;       rho_monte = dblarr(nmonte)
;       phistar_monte = schechter.phistar + randomn(seed,nmonte)*schechter.phistar_err
;       mstar_monte = schechter.mstar + randomn(seed,nmonte)*schechter.mstar_err
;       alpha_monte = schechter.alpha + randomn(seed,nmonte)*schechter.alpha_err
;       for ii = 0, nmonte-1 do begin
;          rho_monte[ii] = total(0.01*10^massaxis*mf_schechter(massaxis,$
;            phistar_monte[ii],mstar_monte[ii],alpha_monte[ii]),/double)
;       endfor
;       schechter.rho_err = djsig(rho_monte)
;    endif    
       
return
end
