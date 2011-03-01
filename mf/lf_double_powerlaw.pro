;+
; NAME:
;   lf_double_powerlaw
; PURPOSE:
;   given double_powerlaw parameters and stellar lumes, return vals
; CAlumlumING SEQUENCE:
;   vals= lf_double_powerlaw(lum, phistar, lstar, alpha)
;     OR
;   vals= lf_double_powerlaw(lum, double_powerlaw)
; INPUTS:
;   lum - [N] set of stellar lumes
;   AND:
;    phistar - phi* parameter in double_powerlaw defn
;    lstar - M* parameter in double_powerlaw defn
;    alpha - alpha parameter ('faint end slope') in double_powerlaw defn
;   OR:
;   double_powerlaw - structure with phistar, lstar, alpha entries
; REVISION HISTORY:
;   14-Apr-2009  Written by John Moustakas, NYU, entirely based
;     on M. Blanton's lf_double_powerlaw
;-

function lf_double_powerlaw, loglum, phistar, loglstar, alpha, beta

    if (n_tags(phistar) gt 0) then begin
       use_phistar = phistar.phistar
       use_loglstar = phistar.lstar
       use_alpha = phistar.alpha
       use_beta = phistar.beta
    endif else begin
       use_phistar = phistar
       use_loglstar = loglstar
       use_alpha = alpha
       use_beta = beta
    endelse

; integrated LF and differential LF (see Yahil+91)
    intlf = use_phistar*(10D^(loglum-use_loglstar))^(-use_alpha)*(1.0+10D^(loglum-use_loglstar))^(-use_beta)
    val = alog(10)*10D^loglum*(use_alpha/10D^loglum+use_beta/(10D^loglum+10D^use_loglstar))*intlf

return, val
end
