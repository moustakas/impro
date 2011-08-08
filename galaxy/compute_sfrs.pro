;+
; NAME:
;       COMPUTE_SFRS()
;
; PURPOSE:
;       Compute emission-line star-formation rates.  This routine must
;       be called after COMPUTE_LINELUMS().  
;
; CALLING SEQUENCE:
;       sfrs = compute_sfrs(linefit)
;
; INPUTS:
;       linefit - can be a structure of either raw or
;                 reddening-corrected fluxes
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       sfrs - 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jan 07, U of A, excised from
;          COMPUTE_LINELUMS() 
;-

function compute_sfrs, linefit

    nspec = n_elements(linefit)
    if (nspec eq 0L) then begin
       print, 'Syntax - sfrs = compute_sfrs(linefit)'
       return, -1L
    endif

    R_HaHb = 2.86 ; intrinsic Ha/Hb ratio
;   R_HaHb = return_tbalmer(/HaHb) ; intrinsic Ha/Hb ratio
    lsunerg = 3.826D33             ; [erg/s]

; compute H-alpha and H-beta SFRs from the reddening-corrected
; luminosities 

    sfrs = {$
      sfr_h_alpha:     -999.0, $
      sfr_h_alpha_err: -999.0, $
      sfr_h_beta:      -999.0, $
      sfr_h_beta_err:  -999.0}
    sfrs = replicate(sfrs,nspec)

    if tag_exist(linefit,'H_ALPHA_LUM') then begin

       good = where(linefit.h_alpha_lum[1] gt 0.0,ngood)
       if (ngood ne 0L) then begin

          lum = 10.0^linefit[good].h_alpha_lum[0]
          lum_err = linefit[good].h_alpha_lum[1]*lum*alog(10.0)
          
          sfr = 7.9D-42*lsunerg*lum ; K98 SFR
          sfr_err = 7.9D-42*lsunerg*lum_err

          sfrs[good].sfr_h_alpha_err = sfr_err/sfr/alog(10.0)
          sfrs[good].sfr_h_alpha = alog10(sfr) 

       endif

    endif

    if tag_exist(linefit,'H_BETA_LUM') then begin
    
       good = where(linefit.h_beta_lum[1] gt 0.0,ngood)
       if (ngood ne 0L) then begin

          lum = R_HaHb*10.0^linefit[good].h_beta_lum[0]
          lum_err = R_HaHb*linefit[good].h_beta_lum[1]*lum*alog(10.0)
          
          sfr = 7.9D-42*lsunerg*lum
          sfr_err = 7.9D-42*lsunerg*lum_err
          
          sfrs[good].sfr_h_beta_err = sfr_err/sfr/alog(10.0)
          sfrs[good].sfr_h_beta = alog10(sfr) 

       endif

    endif

return, sfrs
end    
