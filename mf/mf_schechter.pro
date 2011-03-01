;+
; NAME:
;   mf_schechter
; PURPOSE:
;   given schechter parameters and stellar masses, return vals
; CALLING SEQUENCE:
;   vals= mf_schechter(mass, phistar, mstar, alpha)
;     OR
;   vals= lf_schechter(mass, schechter)
; INPUTS:
;   mass - [N] set of stellar masses
;   AND:
;    phistar - phi* parameter in schechter defn
;    mstar - M* parameter in schechter defn
;    alpha - alpha parameter ('faint end slope') in schechter defn
;   OR:
;   schechter - structure with phistar, mstar, alpha entries
; REVISION HISTORY:
;   14-Apr-2009  Written by John Moustakas, NYU, based largely 
;     on M. Blanton's lf_schechter
;-
;------------------------------------------------------------------------------
function mf_schechter, logmass, phistar, logmstar, alpha
    
; note that MASS and MSTAR and should be logarithms!  (otherwise MPFIT
; in MF_FIT_SCHECHTER has trouble returning the errors)
    
    if (n_tags(phistar) gt 0) then begin
       use_phistar = phistar.phistar
       use_logmstar = phistar.mstar
       use_alpha = phistar.alpha
    endif else begin
       use_phistar = phistar
       use_logmstar = logmstar
       use_alpha = alpha
    endelse

    mratio = 10D^(logmass-use_logmstar)
    val = alog(10)*exp(-mratio) * use_phistar*mratio^use_alpha * mratio

return,val
end
