;+
; NAME:
;   INTEGRATE_MF()
;
; PURPOSE:
;   Integrate a 1/Vmax-weighted stellar mass function and, optionally,
;   the best-fitting Schechter model.
;
; INPUTS: 
;   mf_vmax - 1/Vmax weighted stellar mass function in the style of
;     IM_MF_VMAX()
;
; OPTIONAL INPUTS: 
;   schechter - input MF_SCHECHTER or MF_SCHECHTER_PLUS-style data
;     structure (see /DOUBLE)
;
; KEYWORD PARAMETERS: 
;   double - assume that SCHECHTER describes the double Schechter
;     function used by MF_SCHECHTER_PLUS
;
; OUTPUTS: 
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;   Error checking isn't great.
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 15, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function do_schechter_integral, schechter, double=double

    mass = range(0.0,100.0,500)
    
; rho_tot = Phi**M*Gamma(alpha+2)
    rho = schechter.phistar*10^schechter.logmstar*gamma(schechter.alpha+2)
    
; rho(>M) = Phi**M*IGamma(alpha+2,M/M*)
    rho = schechter.phistar*10^schechter.logmstar*$
      (1-igamma(schechter.alpha+2,10^(9.0-schechter.logmstar)))
    
return, rho
end

function do_mf_integral, mass, phi, phierr=phierr, lomass=lomass, $
  himass=himass, toterr=toterr, nmonte=nmonte

; rho(>M) = 
    rho = im_integral(10^mass,phi/alog(10),10^lomass,10^himass)
    
return, rho
end

function integrate_mf, mf_vmax, schechter=schechter, nmonte=nmonte

    if (n_elements(mf_vmax) eq 0) and (n_elements(schechter) eq 0) then begin
       doc_library, 'integrate_mf'
       return, -1
    endif

    if (n_elements(nmonte) eq 0) then nmonte = 300

; initialize the output data structure
    int = {rho_tot: -1.0, rho_tot_err: -1.0, rho_model_tot: -1.0, rho_model_tot_err: -1.0}


; integrate the data
    if (n_elements(mf_vmax) ne 0) then begin
       gd = where(mf_vmax.limit eq 1,ngd)
       if (ngd ne 0) then int.rho_tot = do_mf_integral(mf_vmax.mass[gd],$
         mf_vmax.phi[gd],phierr=mf_vmax.phi[gd],lomass=min(mf_vmax.mass[gd]),$
         himass=max(mf_vmax.mass[gd]))
    endif

    mrho = do_schechter_integral(schechter,double=double)
    
stop
    
    
; alpha_err is for a very particular application; frequently, the
; faint-end slope alpha is not well-constrained and must be fixed;
; however, when computing the uncertainty on the stellar mass density
; one would like to take the *uncertainty* on alpha into account,
; which can be accomplished using this parameter (in practice, I
; simply overwrite the output error from MPFIT, which is zero is the
; parameter is fixed)

; get the total stellar mass density two ways: by integrating the
; observed MF and also by integrating the *model* from 0-->infty  
    schechter.rho = im_integral(10^logmass,phi)
    rho_monte = dblarr(nmonte)
    for ii = 0, nmonte-1L do begin
       phi_monte = phi + randomn(seed,n_elements(phi))*phierr
       rho_monte[ii] = im_integral(10^logmass,phi_monte)
    endfor
    schechter.rho_err = djsig(rho_monte)

; integrate the model    
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
;      use_covar = covar
       if (schechter.phistar_err gt 0.0) and (schechter.mstar_err gt 0.0) then begin
          rand = mrandomn(seed,use_covar,nmonte)
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
       use_alpha_err = schechter.alpha_err
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
       
    
return, 1
end
    
