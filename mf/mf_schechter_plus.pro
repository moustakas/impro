;+
; NAME:
;   mf_schechter_plus
; PURPOSE:
;   given schechter parameters and stellar masses, return vals
; CALLING SEQUENCE:
;   vals= mf_schechter(mass, phistar, mstar, alpha, phiplus, alphaplus)
;     OR
;   vals= mf_schechter(mass, schechter)
; INPUTS:
;   mass - [N] set of absolute magnitudes
;   AND:
;    phistar - phi* parameter in 1st schechter fn
;    mstar - M* parameter in both schechter fns
;    alpha - alpha parameter ('faint end slope') in 1st schechter fn
;    phiplus - phi* parameter in 2nd schechter defn
;    alphaplus - alpha parameter ('faint end slope') in 2nd schechter fn
;   OR:
;    schechter_plus - structure containing double schechter pars 
;                     (.PHISTAR, .MSTAR, .ALPHA, .PHIPLUS, .ALPHA_PLUS)
; COMMENTS:
;   The double Schechter function is the sum of two Schechter
;   functions with the same MSTAR but different PHISTAR and ALPHA
;   values.
; REVISION HISTORY:
;   20-Oct-2003  Written by Mike Blanton, NYU
;-

function mf_schechter_plus, logmass, phistar, logmstar, $
  alpha, phiplus, alphaplus

; note that MASS and MSTAR and should be logarithms!  (otherwise MPFIT
; in MF_FIT_SCHECHTER has trouble returning the errors)

    if(n_tags(phistar) gt 0) then begin
       use_phistar = phistar.phistar
       use_logmstar = phistar.mstar
       use_alpha = phistar.alpha
       use_phiplus = phistar.phiplus
       use_alphaplus = phistar.alphaplus
    endif else begin
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
       use_phiplus = phiplus
       use_alphaplus = alphaplus
    endelse

    mratio = 10D^(logmass-use_logmstar)
    val = alog(10)*exp(-mratio)*mratio * (use_phistar*mratio^use_alpha + $
      use_phiplus*mratio^use_alphaplus)

return, val
end
