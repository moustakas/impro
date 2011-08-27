;+
; NAME:
;       LM_IGMTAU()
;
; PURPOSE:
;       For a given wavelength and redshift, compute the effective 
;       optical depth due to the IGM, from the Lyman-series opacity
;       and the photoelectric absorption past the Lyman-limit. 
;
; CALLING SEQUENCE:
;       taueff = lm_igmtau(ww,z)
;
; INPUTS:
;       ww - wavelength (Angstroms) [NWW]
;       z  - redshift
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       taueff - effective optical depth [NWW]
;
; MODIFICATION HISTORY:
;       Leonidas A. Moustakas, 2003 March - written
;       J. Moustakas, 2005 Mar 22, U of A, documented, error checking
;          added 
;-

FUNCTION lm_igmtau, ww, z

    nww = n_elements(ww)
    nz = n_elements(z)

    if (nww eq 0L) or (nz eq 0L) then begin
       print, 'Syntax - tauff = lm_igmtau(ww,z)'
       return, 0.0
    endif
    
    z1 = (1.+z)*1.d0
    lymanlimit = 912.d0
    taueff = fltarr(n_elements(ww))
    lyser = [1216.,1026.,973.,950.]*1.d0
    Alyser = [3.6e-3,1.7e-3,1.2e-3,9.3e-4]*1.d0

; series' opacity

    FOR ii=0L, n_elements(lyser)-1L DO $
      taueff = taueff + Alyser[ii]*(ww/lyser[ii])^3.46 * (ww LT lyser[ii]*z1)

; photoelectric absorption past lyman limit

    xc = ww/lymanlimit          ; 1.+zc
    taueff = taueff + $
      (0.25*xc^3*(z1^0.46-xc^0.46)+$
      9.4*xc^1.5*(z1^0.18-xc^0.18)-$
      0.7*xc^3*(xc^(-1.32)-z1^(-1.32))-$
      0.023*(z1^1.68-xc^1.68)) * (ww LE lymanlimit*z1)

; at wavelengths below about 610 Angstroms, the fitting function above
; is not valid.  therefore, force an exponential cutoff...

    wwcut = 650.*z1
    taucut = max(taueff)
    wwcutidx = where(ww LT wwcut)
    IF total(wwcutidx) NE -1 THEN $
      taueff[wwcutidx] = exp(wwcut/ww[wwcutidx]-1)*taucut

;   taueff = taueff + 
;   taueff = taueff + max(taueff)*(ww LT wwcut)
;   taueff = taueff + taueff*exp(-wwcut/ww+1.)*(ww LE wwcut)
    
;   wwbreak = (where(ww EQ lymanlimit*z1))[0]
;   taueff = taueff + (ww LT 700*z1)*(taueff[ww[wwbreak]]*exp(ww/1000.)/exp(ww[wwbreak]))

return, taueff
end
