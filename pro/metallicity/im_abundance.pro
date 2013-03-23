;+
; NAME:
;   IM_ABUNDANCE()
;
; PURPOSE:
;   Compute empirical abundances and various relevant line-ratios
;   from strong emission lines. 
;
; INPUTS:
;   line   - input structure of de-reddened line fluxes
;
; OPTIONAL INPUTS:
;   snrcut_abundance - require that individual emission lines have
;                      a signal-to-noise ratio greater than
;                      SNRCUT_ABUNDANCE 
;   ewalpha          - alpha parameter for Kobulnicky & Phillips
;                      (2003) EW abundances
;   nmonte           - number of Monte Carlo realizations; note
;                      that if NMONTE=0 then the abundance
;                      uncertainties for some calibrations will
;                      not be computed 
;
; KEYWORD PARAMETERS:
;   use_4959  - use the measured 4959 emission-line flux; the
;               default is to assume 4959 = 5007/2.984; note that
;               if LINE has not been reddening-corrected, then the
;               assumption of the intrinsic ratio will be wrong by
;               a little bit 
;   getdensity - compute electron densities (slow!)
;   silent    - do not print messages to STDOUT
;
; OUTPUTS:
;   abund - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Feb 8, U of A, written
;   jm04jul06uofa - added empirical Te-based abundances 
;   jm04jul22uofa - added SILENT and COMBINATION keywords 
;   jm04jul27uofa - added a warning flag if R23>1.1; added
;                   ELECTRONDENSITY keyword
;   jm04nov03uofa - documentation updated; HBETA keyword replaced
;                   with HbHg keyword for consistency with
;                   IUNRED_LINEDUST() 
;   jm04dec08uofa - general development of various empirical
;                   diagnostics; remove objects with log R23 > 1 
;   jm05jan01uofa - compute EW line-ratios and abundances
;   jm05jan06uofa - do not require reddening-corrected line
;                   fluxes
;   jm05feb08uofa - added NOTSTRICT keyword
;   jm05mar31uofa - cleaned up the structure field names
;   jm05sep05uofa - removed NOTSTRICT keyword in favor of
;                   R23STRICT keyword
;   jm06mar10uofa - added EWALPHA and NMONTE parameters 
;   jm07dec10nyu  - added S23
;   jm08feb06nyu  - somewhat major overhaul; cleaned up and
;                   streamlined; removed NORMALIZE and
;                   ELECTRONDENSITY keywords; added JUSTEW and
;                   JUSTFLUX keywords; KK04, M91, and PT05
;                   abundances now computed formally using their
;                   own routines
;   jm08mar20nyu - added NODENSITY and USE_4959 keywords; allow
;                  NMONTE=0 
;   jm08aug21nyu - compute the T04 calibration using EWs
;   jm08oct01nyu - added JUSTRATIOS keyword
;   jm10jul23ucsd - NODENSITY keyword changed to GETDENSITY; CL01 and
;     KD02/combined abundances removed since I never use them
;
; Copyright (C) 2004-2008, 2010, John Moustakas
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

function get_r23lines, line, ew=ew
; support routine for when we call, e.g., MONTE_LOG12OH_KK04
    r23lines = replicate({oii: [-999.0,-999.0], oiii: [-999.0,-999.0], $
      hbeta: [-999.0,-999.0]},n_elements(line))
    if keyword_set(ew) then begin
       r23lines.oii = line.oii_3727_ew
       r23lines.oiii = line.oiii_4959_5007_ew
       r23lines.hbeta = line.h_beta_ew
    endif else begin
       r23lines.oii = line.oii_3727
       r23lines.oiii = line.oiii_4959_5007
       r23lines.hbeta = line.h_beta
    endelse
return, r23lines
end    

function im_abundance, line, snrcut_abundance=snrcut_abundance, $
  ewalpha=ewalpha1, justew=justew, justflux=justflux, justratios=justratios, $
  nmonte=nmonte, getdensity=getdensity, use_4959=use_4959, silent=silent
    
    nspec = n_elements(line)
    if (nspec eq 0L) then begin
       doc_library, 'im_abundance'
       return, -1L
    endif

    oratio = 2.984 ; intrinsic 5007/4959 ratio
    ocor = 1.0 + 1.0/oratio
    
; default: compute abundances based on both EWS and FLUXES 
    
    if (not keyword_set(justew)) and (not keyword_set(justflux)) then begin
       justew = 1
       justflux = 1
    endif
    
    if (n_elements(nmonte) eq 0) then nmonte = 250
    if (n_elements(snrcut_abundance) eq 0) then snrcut_abundance = 5.0
    if (not keyword_set(silent)) then splog, 'S/N > '+$
      string(snrcut_abundance,format='(G0.0)')+', '+$
      'NMONTE = '+string(nmonte,format='(G0.0)')

; initialize the output data structure    
    
    abund = {$
; flux-ratios of interest
      zstrong_r23:                       -999.0, $ ; ([O II] + [O III] 4959,5007) / H-beta
      zstrong_r23_err:                   -999.0, $
      zstrong_o32:                       -999.0, $ ; [O III] 4959,5007 / [O II] 3727
      zstrong_o32_err:                   -999.0, $
      zstrong_p:                         -999.0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      zstrong_p_err:                     -999.0, $
      zstrong_oiiinii:                   -999.0, $ ; ([O III] 5007 / H-beta) / ([N II] 6584 / H-alpha)
      zstrong_oiiinii_err:               -999.0, $
      zstrong_niiha:                     -999.0, $ ; [N II] 6584 / H-alpha
      zstrong_niiha_err:                 -999.0, $
      zstrong_niioii:                    -999.0, $ ; [N II] 6584 / [O II] 3727
      zstrong_niioii_err:                -999.0, $
      zstrong_niisii:                    -999.0, $ ; [N II] 6584 / [S II] 6716,6731
      zstrong_niisii_err:                -999.0, $
      zstrong_oiioiii:                   -999.0, $ ; [O II] 3727 / [O III] 5007
      zstrong_oiioiii_err:               -999.0, $
      zstrong_oiiihb:                    -999.0, $
      zstrong_oiiihb_err:                -999.0, $
      zstrong_s23:                       -999.0, $ ; ([S II] 6716,6731 + [S III] 9069,9532) / H-beta
      zstrong_s23_err:                   -999.0, $

; density-sensitive ratios and densities
      zstrong_sii:                       -999.0, $ ; [S II] 6716 / [S II] 6731
      zstrong_sii_err:                   -999.0, $ 
      zstrong_oii:                       -999.0, $ ; [O II] 3726 / [O II] 3729
      zstrong_oii_err:                   -999.0, $ 
      zstrong_sii_dens:                  -999.0, $ ; [S II] electron density (default 100 cm-3)
      zstrong_sii_dens_err:              -999.0, $
      zstrong_oii_dens:                  -999.0, $ ; [O II] electron density (default 100 cm-3)
      zstrong_oii_dens_err:              -999.0, $

; R23-based O/H calibrations
;     zstrong_12oh_m91_frac:             -999.0, $ ; McGaugh (1991)
      zstrong_12oh_m91_avg:              -999.0, $
      zstrong_12oh_m91_avg_err:          -999.0, $
      zstrong_12oh_m91_upper:            -999.0, $
      zstrong_12oh_m91_upper_err:        -999.0, $
      zstrong_12oh_m91_lower:            -999.0, $
      zstrong_12oh_m91_lower_err:        -999.0, $

;     zstrong_12oh_kk04_frac:            -999.0, $ ; Kobulnicky & Kewley (2004)
      zstrong_12oh_kk04_avg:             -999.0, $
      zstrong_12oh_kk04_avg_err:         -999.0, $
      zstrong_12oh_kk04_upper:           -999.0, $
      zstrong_12oh_kk04_upper_err:       -999.0, $
      zstrong_12oh_kk04_lower:           -999.0, $
      zstrong_12oh_kk04_lower_err:       -999.0, $
      zstrong_logu_kk04_avg:             -999.0, $
      zstrong_logu_kk04_avg_err:         -999.0, $
      zstrong_logu_kk04_upper:           -999.0, $
      zstrong_logu_kk04_upper_err:       -999.0, $
      zstrong_logu_kk04_lower:           -999.0, $
      zstrong_logu_kk04_lower_err:       -999.0, $
      zstrong_converge_kk04_lower:            0, $
      zstrong_converge_kk04_upper:            0, $
                                         
;     zstrong_12oh_pt05_frac:            -999.0, $ ; Pilyugin & Thuan (2005)
      zstrong_12oh_pt05_avg:             -999.0, $
      zstrong_12oh_pt05_avg_err:         -999.0, $
      zstrong_12oh_pt05_upper:           -999.0, $
      zstrong_12oh_pt05_upper_err:       -999.0, $
      zstrong_12oh_pt05_lower:           -999.0, $
      zstrong_12oh_pt05_lower_err:       -999.0, $

      zstrong_12oh_zkh94:                -999.0, $ ; ZKH94 (upper branch only)
      zstrong_12oh_zkh94_err:            -999.0, $
      zstrong_12oh_t04:                  -999.0, $ ; Tremonti et al. 2004 (upper branch only)
      zstrong_12oh_t04_err:              -999.0, $

; other O/H calibrations
      zstrong_12oh_niiha_pettini:        -999.0, $ ; Pettini & Pagel (2004) [N II]/Ha calibration
      zstrong_12oh_niiha_pettini_err:    -999.0, $ 
      zstrong_12oh_oiiinii_pettini:      -999.0, $ ; empirical abundance based on the Pettini & Pagel (2004) calibration
      zstrong_12oh_oiiinii_pettini_err:  -999.0, $ 
;     zstrong_12oh_cl01A:                -999.0, $ ; Charlot & Longhetti (2001), Case A
;     zstrong_12oh_cl01A_err:            -999.0, $
      zstrong_12oh_kd02_niioii:          -999.0, $ ; Kewley & Dopita (2002) using [N II]/[O II]
      zstrong_12oh_kd02_niioii_err:      -999.0, $
;     zstrong_12oh_kd02_combined:        -999.0, $ ; Kewley & Dopita (2002) combined abundance diagnostic
;     zstrong_12oh_kd02_combined_err:    -999.0, $

; EW-ratios of interest
      ewalpha:                              1.0, $
      zstrong_ew_r23:                    -999.0, $ ; {ew([O II]) + ew([O III] 4959,5007)} / ew(H-beta)
      zstrong_ew_r23_err:                -999.0, $
      zstrong_ew_o32:                    -999.0, $ ; ew([O III] 4959,5007) / ew([O II] 3727)
      zstrong_ew_o32_err:                -999.0, $
      zstrong_ew_p:                      -999.0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      zstrong_ew_p_err:                  -999.0, $

; EW(R23)-based O/H calibrations
;     zstrong_ew_12oh_m91_frac:          -999.0, $ ; McGaugh (1991)
      zstrong_ew_12oh_m91_avg:           -999.0, $
      zstrong_ew_12oh_m91_avg_err:       -999.0, $
      zstrong_ew_12oh_m91_upper:         -999.0, $
      zstrong_ew_12oh_m91_upper_err:     -999.0, $
      zstrong_ew_12oh_m91_lower:         -999.0, $
      zstrong_ew_12oh_m91_lower_err:     -999.0, $
                                         
;     zstrong_ew_12oh_kk04_frac:         -999.0, $ ; Kobulnicky & Kewley (2004)
      zstrong_ew_12oh_kk04_avg:          -999.0, $
      zstrong_ew_12oh_kk04_avg_err:      -999.0, $
      zstrong_ew_12oh_kk04_upper:        -999.0, $
      zstrong_ew_12oh_kk04_upper_err:    -999.0, $
      zstrong_ew_12oh_kk04_lower:        -999.0, $
      zstrong_ew_12oh_kk04_lower_err:    -999.0, $
      zstrong_ew_logu_kk04_avg:          -999.0, $
      zstrong_ew_logu_kk04_avg_err:      -999.0, $
      zstrong_ew_logu_kk04_upper:        -999.0, $
      zstrong_ew_logu_kk04_upper_err:    -999.0, $
      zstrong_ew_logu_kk04_lower:        -999.0, $
      zstrong_ew_logu_kk04_lower_err:    -999.0, $
      zstrong_ew_converge_kk04_lower:         0, $
      zstrong_ew_converge_kk04_upper:         0, $
                                         
;;     zstrong_ew_12oh_pt05_frac:         -999.0, $ ; Pilyugin & Thuan (2005)
;      zstrong_ew_12oh_pt05_avg:          -999.0, $
;      zstrong_ew_12oh_pt05_avg_err:      -999.0, $
;      zstrong_ew_12oh_pt05_upper:        -999.0, $
;      zstrong_ew_12oh_pt05_upper_err:    -999.0, $
;      zstrong_ew_12oh_pt05_lower:        -999.0, $
;      zstrong_ew_12oh_pt05_lower_err:    -999.0, $
      
      zstrong_ew_12oh_T04:               -999.0, $ ; Tremonti et al. 2004 (upper branch only)
      zstrong_ew_12oh_T04_err:           -999.0}

    abund = replicate(abund,nspec)

; ###########################################################################
; compute EW-ratios and abundances
; ###########################################################################

    if keyword_set(justew) then begin

       if (not keyword_set(silent)) then begin
          if keyword_set(justratios) then splog, 'Computing EW line-ratios' else $
            splog, 'Computing line-ratios and abundances from EWs:'
       endif
       
; choose EWALPHA and then "correct" EW([O II]) to account for it
       
       if (n_elements(ewalpha1) eq 0L) then begin
          splog, 'Using alpha = 1.0 for all galaxies'
          ewalpha = replicate(1.0,nspec) 
       endif else begin
          if (n_elements(ewalpha1) eq 1L) then begin
             splog, 'Using alpha = '+string(ewalpha1,format='(G0.0)')+' for all galaxies'
             ewalpha = replicate(ewalpha1,nspec) 
          endif else begin
             if (nspec ne n_elements(ewalpha1)) then begin
                splog, 'Dimensions of LINE and EWALPHA do not agree'
                return, abund
             endif else begin
                splog, 'Using individual alpha value for each galaxy'
                ewalpha = ewalpha1
             endelse
          endelse
       endelse
       abund.ewalpha = ewalpha
       
       ewline = line
       if (tag_exist(ewline[0],'OII_3727_EW')) then begin
          doit = where((ewline.oii_3727_ew[1] gt 0.0),ndoit)
          if (ndoit ne 0L) then begin
             ewline[doit].oii_3727_ew[0] = ewalpha*ewline[doit].oii_3727_ew[0]
             ewline[doit].oii_3727_ew[1] = ewalpha*ewline[doit].oii_3727_ew[1]
          endif
       endif

; if /USE_4959 then use the measured [O III] 4959 line, otherwise use
; the intrinsic ratio, 2.984; note that if LINE has not been
; reddening-corrected, then the assumption of the intrinsic ratio will
; be wrong by a little bit (although it doesn't really matter for EWs)
       ewline = struct_addtags(temporary(ewline),replicate({oiii_4959_5007_ew: [0.0,-2.0]},nspec))
       if keyword_set(use_4959) then begin
          if tag_exist(ewline[0],'OIII_4959_EW') and tag_exist(ewline[0],'OIII_5007_EW') then begin
             doit = where((ewline.oiii_4959_ew[1] gt 0.0) and (ewline.oiii_5007_ew[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then begin
                ewline[doit].oiii_4959_5007_ew[0] = ewline[doit].oiii_4959_ew[0]+ewline[doit].oiii_5007_ew[0]
                ewline[doit].oiii_4959_5007_ew[1] = sqrt(ewline[doit].oiii_4959_ew[1]^2+ewline[doit].oiii_5007_ew[1]^2)
             endif
          endif 
       endif else begin
          if (tag_exist(ewline[0],'OIII_5007_EW')) then begin
             doit = where((ewline.oiii_5007_ew[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then ewline[doit].oiii_4959_5007_ew = ocor*ewline[doit].oiii_5007_ew
          endif 
       endelse
       
; EW(R23)
       if (tag_exist(ewline[0],'OII_3727_EW') and tag_exist(ewline[0],'OIII_4959_5007_EW') and $
         tag_exist(ewline[0],'H_BETA_EW')) then begin
          lineratio, ewline, ['OII_3727_EW','OIII_4959_5007_EW'], 'H_BETA_EW', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_r23 = x
             abund[index].zstrong_ew_r23_err = xerr
          endif
       endif

; EW(O32)
       if (tag_exist(ewline[0],'OII_3727_EW') and tag_exist(ewline[0],'OIII_4959_5007_EW')) then begin
          lineratio, ewline, 'OIII_4959_5007_EW', 'OII_3727_EW', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_o32 = x
             abund[index].zstrong_ew_o32_err = xerr
          endif
       endif

; EW(P)-parameter (Pilyugin 2001)
       if tag_exist(ewline[0],'OIII_4959_5007_EW') and tag_exist(ewline[0],'OII_3727_EW') then begin
          lineratio, ewline, 'OIII_4959_5007_EW', ['OII_3727_EW','OIII_4959_5007_EW'], $
            '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_ew_p = x
             abund[index].zstrong_ew_p_err = xerr
          endif
       endif

; save memory!
       ewr23lines = get_r23lines(ewline,/ew)
       ewline = 0

       if (keyword_set(justratios) eq 0) then begin ; compute abundances
          
; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999) - based on EWs
          if (not keyword_set(silent)) then splog, '   M91...'
          index = where((abund.zstrong_ew_o32 gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_m91(abund[index].zstrong_ew_r23,$
               abund[index].zstrong_ew_o32,line=ewr23lines[index],$
               nmonte=nmonte)

;            abund[index].zstrong_ew_12oh_m91_frac = result.frac
             abund[index].zstrong_ew_12oh_m91_avg    = result.log12oh_avg
             abund[index].zstrong_ew_12oh_m91_upper  = result.log12oh_upper
             abund[index].zstrong_ew_12oh_m91_lower  = result.log12oh_lower

             abund[index].zstrong_ew_12oh_m91_avg_err    = result.log12oh_avg_err
             abund[index].zstrong_ew_12oh_m91_upper_err  = result.log12oh_upper_err
             abund[index].zstrong_ew_12oh_m91_lower_err  = result.log12oh_lower_err

          endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - EW/R23
          if (not keyword_set(silent)) then splog, '   KK04...'
          index = where((abund.zstrong_ew_o32 gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_kk04(abund[index].zstrong_ew_r23,$
               abund[index].zstrong_ew_o32,line=ewr23lines[index],$
               nmonte=nmonte)

;            abund[index].zstrong_ew_12oh_kk04_frac = result.frac
             abund[index].zstrong_ew_12oh_kk04_avg   = result.log12oh_avg
             abund[index].zstrong_ew_12oh_kk04_upper = result.log12oh_upper
             abund[index].zstrong_ew_12oh_kk04_lower = result.log12oh_lower
             abund[index].zstrong_ew_logu_kk04_avg   = result.logu_avg
             abund[index].zstrong_ew_logu_kk04_upper = result.logu_upper
             abund[index].zstrong_ew_logu_kk04_lower = result.logu_lower

             abund[index].zstrong_ew_converge_kk04_upper = result.converge_upper
             abund[index].zstrong_ew_converge_kk04_lower = result.converge_lower
             
             abund[index].zstrong_ew_12oh_kk04_avg_err   = result.log12oh_avg_err
             abund[index].zstrong_ew_12oh_kk04_upper_err = result.log12oh_upper_err
             abund[index].zstrong_ew_12oh_kk04_lower_err = result.log12oh_lower_err
             abund[index].zstrong_ew_logu_kk04_avg_err   = result.logu_avg_err
             abund[index].zstrong_ew_logu_kk04_upper_err = result.logu_upper_err
             abund[index].zstrong_ew_logu_kk04_lower_err = result.logu_lower_err
          endif

;; ---------------------------------------------------------------------------    
;; Pilyugin & Thuan (2005) - EWs
;          if (not keyword_set(silent)) then splog, '   PT05...'
;          index = where((abund.zstrong_ew_p gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
;          if (nindex ne 0L) then begin
;
;             result = monte_log12oh_pt05(abund[index].zstrong_ew_r23,$
;               abund[index].zstrong_ew_p,line=ewr23lines[index],nmonte=nmonte)
;
;;            abund[index].zstrong_ew_12oh_pt05_frac = result.frac
;             abund[index].zstrong_ew_12oh_pt05_avg    = result.log12oh_avg
;             abund[index].zstrong_ew_12oh_pt05_upper  = result.log12oh_upper
;             abund[index].zstrong_ew_12oh_pt05_lower  = result.log12oh_lower
;
;             abund[index].zstrong_ew_12oh_pt05_avg_err    = result.log12oh_avg_err
;             abund[index].zstrong_ew_12oh_pt05_upper_err  = result.log12oh_upper_err
;             abund[index].zstrong_ew_12oh_pt05_lower_err  = result.log12oh_lower_err
;          endif
          
; ---------------------------------------------------------------------------    
; Tremonti et al. (2004)
          if (keyword_set(silent) eq 0) then splog, '   T04...'
          index = where(abund.zstrong_ew_r23 gt -900.0,nindex)
          if (nindex ne 0L) then begin
             result = monte_log12oh_t04(abund[index].zstrong_ew_r23,$
               line=ewr23lines[index],nmonte=nmonte)
             abund[index].zstrong_ew_12oh_t04     = result.log12oh
             abund[index].zstrong_ew_12oh_t04_err = result.log12oh_err
          endif
       endif                    ; close the JUSTRATIOS case
    endif                       ; close the JUSTEW case

; ###########################################################################
; compute flux-ratios and abundances
    if keyword_set(justflux) then begin
       if (keyword_set(silent) eq 0) then begin
          if keyword_set(justratios) then splog, 'Computing line-flux ratios' else $
            splog, 'Computing line-ratios and abundances from fluxes:'
       endif

; if /USE_4959 then use the measured [O III] 4959 line, otherwise use
; the intrinsic ratio, 2.984; note that if LINE has not been
; reddening-corrected, then the assumption of the intrinsic ratio will
; be wrong by a little bit

       fluxline = struct_addtags(line,replicate({oiii_4959_5007: [0.0,-2.0]},nspec))
       if keyword_set(use_4959) then begin
          if tag_exist(fluxline[0],'OIII_4959') and tag_exist(fluxline[0],'OIII_5007') then begin
             doit = where((fluxline.oiii_4959[1] gt 0.0) and (fluxline.oiii_5007[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then begin
                fluxline[doit].oiii_4959_5007[0] = fluxline[doit].oiii_4959[0]+fluxline[doit].oiii_5007[0]
                fluxline[doit].oiii_4959_5007[1] = sqrt(fluxline[doit].oiii_4959[1]^2+fluxline[doit].oiii_5007[1]^2)
             endif
          endif 
       endif else begin
          if (tag_exist(fluxline[0],'OIII_5007')) then begin
             doit = where((fluxline.oiii_5007[1] gt 0.0),ndoit)
             if (ndoit ne 0L) then fluxline[doit].oiii_4959_5007 = ocor*fluxline[doit].oiii_5007
          endif 
       endelse

; R23
       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_4959_5007') and $
         tag_exist(fluxline[0],'H_BETA')) then begin
          lineratio, fluxline, ['OII_3727','OIII_4959_5007'], 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_r23 = x
             abund[index].zstrong_r23_err = xerr
          endif
       endif 

; O32
       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_4959_5007')) then begin
          lineratio, fluxline, 'OIII_4959_5007', 'OII_3727', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_o32 = x
             abund[index].zstrong_o32_err = xerr
          endif
       endif 

; P-parameter (Pilyugin 2001)
       if tag_exist(fluxline[0],'OIII_4959_5007') and tag_exist(fluxline[0],'OII_3727') then begin
          lineratio, fluxline, 'OIII_4959_5007', ['OII_3727','OIII_4959_5007'], $
            '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_p = x
             abund[index].zstrong_p_err = xerr
          endif
       endif 

; ([O III]/Hb)/([N II]/Ha) - Pettini & Pagel (2004)
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'OIII_5007') and $
         tag_exist(fluxline[0],'H_BETA') and tag_exist(fluxline[0],'H_ALPHA')) then begin
          lineratio, fluxline, 'OIII_5007', 'H_BETA', 'NII_6584', 'H_ALPHA', $
            x1, x1err, x2, x2err, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_oiiinii = x1 - x2
             abund[index].zstrong_oiiinii_err = sqrt(x1err^2 + x2err^2)
          endif
       endif

; [N II]/Ha    
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'H_ALPHA')) then begin
          lineratio, fluxline, 'NII_6584', 'H_ALPHA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_niiha = x
             abund[index].zstrong_niiha_err = xerr
          endif
       endif

; [N II]/[O II]
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'OII_3727')) then begin
          lineratio, fluxline, 'NII_6584', 'OII_3727', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_niioii = x
             abund[index].zstrong_niioii_err = xerr
          endif
       endif

; [N II]/[S II]
       if (tag_exist(fluxline[0],'NII_6584') and tag_exist(fluxline[0],'SII_6716') and $
         tag_exist(fluxline[0],'SII_6731')) then begin
          lineratio, fluxline, 'NII_6584', ['SII_6716','SII_6731'], '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_niisii = x
             abund[index].zstrong_niisii_err = xerr
          endif
       endif

; [O II]/[O III]
       if (tag_exist(fluxline[0],'OII_3727') and tag_exist(fluxline[0],'OIII_5007')) then begin
          lineratio, fluxline, 'OII_3727', 'OIII_5007', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_oiioiii = x
             abund[index].zstrong_oiioiii_err = xerr
          endif
       endif

; [O III]/H-beta
       if (tag_exist(fluxline[0],'OIII_5007') and tag_exist(fluxline[0],'H_BETA')) then begin
          lineratio, fluxline, 'OIII_5007', 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
          if (nindex ne 0L) then begin
             abund[index].zstrong_oiiihb = x
             abund[index].zstrong_oiiihb_err = xerr
          endif
       endif

; S23
       if (tag_exist(fluxline[0],'SII_6716') and tag_exist(fluxline[0],'SII_6731') and $
         tag_exist(fluxline[0],'SIII_9069') and tag_exist(fluxline[0],'SIII_9532') and $
         tag_exist(fluxline[0],'H_BETA')) then begin
          lineratio, fluxline, ['SII_6716','SII_6731','SIII_9069','SIII_9532'], 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_S23 = x
             abund[index].zstrong_S23_err = xerr
          endif
       endif

; [S II] 6716 / [S II] 6731
       if (tag_exist(fluxline[0],'SII_6716') and tag_exist(fluxline[0],'SII_6731')) then begin
          lineratio, fluxline, 'SII_6716', 'SII_6731', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_sii = x
             abund[index].zstrong_sii_err = xerr
          endif
       endif

; [O II] 3726 / [O II] 3729
       if (tag_exist(fluxline[0],'OII_3726') and tag_exist(fluxline[0],'OII_3729')) then begin
          lineratio, fluxline, 'OII_3726', 'OII_3729', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
          if (nindex ne 0L) then begin
             abund[index].zstrong_oii = x
             abund[index].zstrong_oii_err = xerr
          endif
       endif

       r23lines = get_r23lines(fluxline)
       fluxline = 0 ; save memory
       
       if (keyword_set(justratios) eq 0) then begin ; compute abundances

; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999)
          if (not keyword_set(silent)) then splog, '   M91...'
          index = where((abund.zstrong_o32 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_m91(abund[index].zstrong_r23,$
               abund[index].zstrong_o32,line=r23lines[index],$
               nmonte=nmonte)

;            abund[index].zstrong_12oh_m91_frac = result.frac
             abund[index].zstrong_12oh_m91_avg    = result.log12oh_avg
             abund[index].zstrong_12oh_m91_upper  = result.log12oh_upper
             abund[index].zstrong_12oh_m91_lower  = result.log12oh_lower

             abund[index].zstrong_12oh_m91_avg_err    = result.log12oh_avg_err
             abund[index].zstrong_12oh_m91_upper_err  = result.log12oh_upper_err
             abund[index].zstrong_12oh_m91_lower_err  = result.log12oh_lower_err

          endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - Flux/R23
          if (keyword_set(silent) eq 0) then splog, '   KK04...'
          index = where((abund.zstrong_o32 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_kk04(abund[index].zstrong_r23,$
               abund[index].zstrong_o32,line=r23lines[index],$
               nmonte=nmonte)
             
;            abund[index].zstrong_12oh_kk04_frac = result.frac
             abund[index].zstrong_12oh_kk04_avg   = result.log12oh_avg
             abund[index].zstrong_12oh_kk04_upper = result.log12oh_upper
             abund[index].zstrong_12oh_kk04_lower = result.log12oh_lower
             abund[index].zstrong_logu_kk04_avg   = result.logu_avg
             abund[index].zstrong_logu_kk04_upper = result.logu_upper
             abund[index].zstrong_logu_kk04_lower = result.logu_lower
             
             abund[index].zstrong_converge_kk04_upper = result.converge_upper
             abund[index].zstrong_converge_kk04_lower = result.converge_lower

             abund[index].zstrong_12oh_kk04_avg_err   = result.log12oh_avg_err
             abund[index].zstrong_12oh_kk04_upper_err = result.log12oh_upper_err
             abund[index].zstrong_12oh_kk04_lower_err = result.log12oh_lower_err
             abund[index].zstrong_logu_kk04_avg_err   = result.logu_avg_err
             abund[index].zstrong_logu_kk04_upper_err = result.logu_upper_err
             abund[index].zstrong_logu_kk04_lower_err = result.logu_lower_err
          endif

; ---------------------------------------------------------------------------    
; Pilyugin & Thuan (2005)
          if (not keyword_set(silent)) then splog, '   PT05...'
          index = where((abund.zstrong_o32 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin

             result = monte_log12oh_pt05(abund[index].zstrong_r23,$
               abund[index].zstrong_p,line=r23lines[index],nmonte=nmonte)

;            abund[index].zstrong_12oh_pt05_frac = result.frac
             abund[index].zstrong_12oh_pt05_avg    = result.log12oh_avg
             abund[index].zstrong_12oh_pt05_upper  = result.log12oh_upper
             abund[index].zstrong_12oh_pt05_lower  = result.log12oh_lower

             abund[index].zstrong_12oh_pt05_avg_err    = result.log12oh_avg_err
             abund[index].zstrong_12oh_pt05_upper_err  = result.log12oh_upper_err
             abund[index].zstrong_12oh_pt05_lower_err  = result.log12oh_lower_err
          endif

; ---------------------------------------------------------------------------    
; Zaritsky et al. (1994)
          index = where((abund.zstrong_r23 gt -900.0),nindex)
          if (nindex ne 0L) then begin
             result = monte_log12oh_zkh94(abund[index].zstrong_r23,$
               line=r23lines[index],nmonte=nmonte)
             abund[index].zstrong_12oh_zkh94     = result.log12oh
             abund[index].zstrong_12oh_zkh94_err = result.log12oh_err
          endif

; ---------------------------------------------------------------------------    
; Tremonti et al. (2004)
          index = where(abund.zstrong_r23 gt -900.0,nindex)
          if (nindex ne 0L) then begin
             result = monte_log12oh_t04(abund[index].zstrong_r23,$
               line=r23lines[index],nmonte=nmonte)
             abund[index].zstrong_12oh_t04     = result.log12oh
             abund[index].zstrong_12oh_t04_err = result.log12oh_err
          endif
          
; ---------------------------------------------------------------------------    
; Pettini & Pagel (2004) [N II]/Ha calibration
          if (not keyword_set(silent)) then splog, '   N2...'
          index = where((abund.zstrong_niiha ge -2.5) and (abund.zstrong_niiha lt -0.3),nindex)
          if (nindex ne 0L) then begin

             c = [8.90D,0.57D]
;            cerr = [0.18D,0.0D]
             cerr = [0.0D,0.0D]
             
             x = abund[index].zstrong_niiha
             xerr = abund[index].zstrong_niiha_err
             
             abund[index].zstrong_12oh_niiha_pettini = c[0] + c[1] * x
             abund[index].zstrong_12oh_niiha_pettini_err = sqrt( cerr[0]^2 + $
               (c[1] * xerr)^2.0 + (cerr[1] * x)^2.0 )

             x = 0 & xerr = 0   ; save memory

          endif
          
; ---------------------------------------------------------------------------    
; Pettini & Pagel (2004) ([O III]/Hb)/([N II]/Ha) calibration
          if (not keyword_set(silent)) then splog, '   O3N2...'
          index = where((abund.zstrong_oiiinii ge -1.0) and (abund.zstrong_oiiinii le 1.9),nindex)
          if (nindex ne 0L) then begin

             c = [8.73D,-0.32D]
             cerr = [0.0D,0.0D]
;         cerr = [0.14D,0.0D]

             x = abund[index].zstrong_oiiinii
             xerr = abund[index].zstrong_oiiinii_err

             abund[index].zstrong_12oh_oiiinii_pettini = c[0] + c[1] * x
             abund[index].zstrong_12oh_oiiinii_pettini_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

             x = 0 & xerr = 0   ; save memory

          endif

;; ---------------------------------------------------------------------------    
;; Charlot & Longhetti (2001, Case A)    
;          if (not keyword_set(silent)) then splog, '   CL01/Case A...'
;          index = where((abund.zstrong_niisii gt -900.0) and (abund.zstrong_oiioiii gt -900.0),nindex)
;          if (nindex ne 0L) then begin
;
;             logx2 = abund[index].zstrong_oiioiii
;             logx2_err = abund[index].zstrong_oiioiii_err
;
;             x2 = 10^logx2/1.5
;             x2_err = alog(10.0)*logx2_err*x2
;
;             logx6 = abund[index].zstrong_niisii
;             logx6_err = abund[index].zstrong_niisii_err
;
;             x6 = 10^logx6/0.85
;             x6_err = alog(10.0)*logx6_err*x6
;
;             term = 5.09D-4 * x2^0.17 * x6^1.17
;             term_err = 5.09D-4 * sqrt( (0.17*x2^(0.17-1)*x6^1.17*x2_err)^2 + (1.17*x6^(1.17-1)*x2^0.17*x6_err)^2 )
;             
;             abund[index].zstrong_12oh_cl01A = 12.0 + alog10(term)
;             abund[index].zstrong_12oh_cl01A_err = term_err/term/alog(10.0)
;             
;             logx2 = 0 & logx2_err = 0 & x2 = 0 & x2_err = 0 & logx6 = 0 & logx6_err = 0 ; save memory
;             x6 = 0 & x6_err = 0 & term = 0 & term_err = 0             
;          endif
          
; ---------------------------------------------------------------------------    
; Kewley & Dopita (2002) [N II]/[O II], use the coefficients from
; Table 3 at the mid-point of the ionization parameter range, log U =
; -3.17 (q=2E7), which should give the same abundance as equation 5
; above 8.6. 
          if (not keyword_set(silent)) then splog, '   KD02/N2O2...'
          index = where((abund.zstrong_niioii gt -900.0),nindex)
          if (nindex ne 0L) then begin

             x = abund[index].zstrong_niioii
             xerr = abund[index].zstrong_niioii_err

             if (nmonte gt 0L) then begin
                xmonte = rebin(reform(x,1,nindex),nmonte,nindex)+$
                  randomn(seed,nmonte,nindex)*rebin(reform(xerr,1,nindex),nmonte,nindex)
             endif

             basecoeff = [1106.87D,-532.154D,96.3733D,-7.81061D,0.239282D]
             coeff = basecoeff # (dblarr(nindex)+1) ; q = 2E7
             coeff[0,*] = coeff[0,*] - x

; find the roots of this polynomial; the correct one must be real and
; between 7.5 and 9.4

             for iobj = 0L, nindex-1L do begin
                roots = fz_roots(coeff[*,iobj],/double)
                for jiter = 0L, 3L do begin ; iterate
                   if (imaginary(roots[jiter]) eq 0.0) then begin
                      if ((abs(roots[jiter]) ge 7.5) and (abs(roots[jiter]) le 9.4)) then $
                        abund[index[iobj]].zstrong_12oh_kd02_niioii = abs(roots[jiter])
                   endif
                endfor
; compute the error
                if (nmonte gt 0L) and (abund[index[iobj]].zstrong_12oh_kd02_niioii gt -900.0) then begin
                   montecoeff = basecoeff # (dblarr(nmonte)+1) ; q = 2E7
                   montecoeff[0,*] = montecoeff[0,*] - xmonte[*,iobj]
                   ohmonte = fltarr(nmonte)-999.0
                   for imonte = 0L, nmonte-1L do begin
                      monteroots = fz_roots(montecoeff[*,imonte],/double)
                      for jiter = 0L, 3L do begin ; iterate
                         if (imaginary(monteroots[jiter]) eq 0.0) then begin
                            if ((abs(monteroots[jiter]) ge 7.5) and (abs(monteroots[jiter]) le 9.4)) then $
                              ohmonte[imonte] = abs(monteroots[jiter])
                         endif
                      endfor
                   endfor
                   goodmonte = where((ohmonte gt -900.0),ngoodmonte)
                   if (ngoodmonte ne 0L) then $
                     abund[index[iobj]].zstrong_12oh_kd02_niioii_err = djsig(ohmonte[goodmonte])
                endif
             endfor 
             
             x = 0 & xerr = 0 & xmonte = 0
             
          endif

;; ---------------------------------------------------------------------------    
;; Kewley & Dopita (2002) combined diagnostic; use KD02 [N II]/[O II]
;; if [N II]/[O II] gives 8.6 < 12 + log (O/H); use the average of M91
;; and ZKH94 if 8.5 < 12 + log (O/H) < 8.6; finally, use the average of
;; CL01A and KK04-R23 if log(O/H)+12 < 8.5.
;          if (not keyword_set(silent)) then splog, '   KD02/Combined...'
;          index = where((abund.zstrong_12oh_kd02_niioii gt -900.0) and (abund.zstrong_12oh_m91_upper gt -900.0) and $
;            (abund.zstrong_12oh_zkh94 gt -900.0) and (abund.zstrong_12oh_kk04_lower gt -900.0) and $
;            (abund.zstrong_12oh_cl01a gt -900.0),nindex)
;
;          if (nindex ne 0L) then begin
;
;             m91Z94_average = (abund[index].zstrong_12oh_m91_upper + abund[index].zstrong_12oh_zkh94) / 2.0
;             kd02C01_average = (abund[index].zstrong_12oh_kk04_lower + abund[index].zstrong_12oh_cl01a) / 2.0
;
;             lo = where(abund[index].zstrong_12oh_kd02_niioii le 8.6,nlo,comp=up,ncomp=nup)
;
;             if (nlo ne 0L) then begin
;
;                upper = where(m91Z94_average[lo] ge 8.5,nupper,comp=lower,ncomp=nlower)
;                if (nupper ne 0L) then abund[index[lo[upper]]].zstrong_12oh_kd02_combined = m91Z94_average[lo[upper]]
;                if (nlower ne 0L) then abund[index[lo[lower]]].zstrong_12oh_kd02_combined = kd02C01_average[lo[lower]]
;                
;             endif
;
;             if (nup ne 0L) then begin
;                abund[index[up]].zstrong_12oh_kd02_combined = abund[index[up]].zstrong_12oh_kd02_niioii
;             endif             
;          endif

; [S II] electron density [cm-3]
          index = where((abund.zstrong_sii gt -900.0),nindex)
          if (nindex ne 0L) and keyword_set(getdensity) then begin
             if (not keyword_set(silent)) then splog, 'Computing [S II] electron densities...'
             temden = im_temden('s_ii',abund[index].zstrong_sii,$
               ratio_err=abund[index].zstrong_sii_err,nmonte=nmonte)
             abund[index].zstrong_sii_dens     = temden.dens
             abund[index].zstrong_sii_dens_err = temden.dens_err
          endif

; [O II] electron density [cm-3]
          index = where((abund.zstrong_oii gt -900.0),nindex)
          if (nindex ne 0L) and keyword_set(getdensity) then begin
             if (not keyword_set(silent)) then splog, 'Computing [O II] electron densities...'
             temden = im_temden('o_ii_dens',abund[index].zstrong_oii,$
               ratio_err=abund[index].zstrong_oii_err,nmonte=nmonte)
             abund[index].zstrong_oii_dens     = temden.dens
             abund[index].zstrong_oii_dens_err = temden.dens_err
          endif
       endif                    ; close the JUSTRATIOS case
    endif                       ; close the JUSTFLUX case
       
return, abund
end

;;;     zstrong_ew_niiha:                  -999.0, $ ; ew([N II] 6584) / ew(H-alpha)
;;;     zstrong_ew_niiha_err:              -999.0, $
;;;     zstrong_oiiihbeta:                 -999.0, $ ; [O III] 5007 / H-beta
;;;     zstrong_oiiihbeta_err:             -999.0, $
;;;     zstrong_r3:                        -999.0, $ ; [O III] 4959,5007 / H-beta
;;;     zstrong_r3_err:                    -999.0, $
;;;     zstrong_T_oiii:                    -999.0, $ ; empirical T(O++) based on [N II]/Ha
;;;     zstrong_T_oiii_err:                -999.0, $
;;;     zstrong_T_oii:                     -999.0, $ ; empirical T(O+) from T(O++) and photoionization models
;;;     zstrong_T_oii_err:                 -999.0, $
;;;     zstrong_log12oh:                   -999.0, $ ; empirical 12+log (O/H) based on empirical T(O++)
;;;     zstrong_log12oh_err:               -999.0, $
;;;     zstrong_log_no:                    -999.0, $ ; empirical log (N/O) based on empirical T(O++)
;;;     zstrong_log_no_err:                -999.0, $
;;;     zstrong_12oh_niiha_denicolo:       -999.0, $ ; Denicolo (2002) [N II]/Ha calibration
;;;     zstrong_12oh_niiha_denicolo_err:   -999.0, $ 
;;;     zstrong_12oh_niiha_moustakas:      -999.0, $ ; empirical abundance based my [N II]/Ha - Te relation
;;;     zstrong_12oh_niiha_moustakas_err:  -999.0, $ 
;;;     zstrong_12oh_niiha_melbourne:      -999.0, $ ; Melbourne & Salzer (2002) [N II]/Ha calibration
;;;     zstrong_12oh_niiha_melbourne_err:  -999.0, $ 
;;;     zstrong_12oh_o32:                  -999.0, $ ; empirical abundance based my o32 - Te relation
;;;     zstrong_12oh_o32_err:              -999.0, $ 
;;;     zstrong_12oh_P:                    -999.0, $ ; empirical abundance based my P - Te relation
;;;     zstrong_12oh_P_err:                -999.0, $ 
;;;     zstrong_12oh_oiiinii:              -999.0, $ ; empirical abundance based my ([O III]/Hb)/([N II]/Ha) calibration - WRONG!!!
;;;     zstrong_12oh_oiiinii_err:          -999.0, $ 
;;;     zstrong_12oh_oiiinii_moustakas:    -999.0, $ ; empirical abundance based on the Pettini & Pagel (2004) calibration
;;;     zstrong_12oh_oiiinii_moustakas_err:-999.0, $ 
;;;     zstrong_12oh_cl01F:                -999.0, $ ; Charlot & Longhetti (2001), Case F
;;;     zstrong_12oh_cl01F_err:            -999.0, $
;;;     zstrong_12oh_m91_contini:          -999.0, $ ; Contini et al. (2002) procedure for m91
;;;     zstrong_12oh_m91_contini_err:      -999.0, $
;;;     zstrong_12oh_m91_kd02:             -999.0, $ ; Kewley & Dopita (2002) procedure for m91
;;;     zstrong_12oh_m91_kd02_err:         -999.0, $
;;;     zstrong_12oh_m91_o32:              -999.0, $ ; use o32 to choose the branch
;;;     zstrong_12oh_m91_o32_err:          -999.0, $
;;;     zstrong_12oh_kd02_r23:             -999.0, $ ; Kewley & Dopita (2002) using r23
;;;     zstrong_12oh_kd02_r23_err:         -999.0, $
;;;     zstrong_12oh_pt05_o32:             -999.0, $ ; use o32 to choose the branch
;;;     zstrong_12oh_pt05_o32_err:         -999.0, $
;;;     zstrong_12oh_pt05_average:         -999.0, $ ; average of the lower and upper calibrations, to be used in the turn-around region
;;;     zstrong_12oh_pt05_average_err:     -999.0, $
;;;     zstrong_12oh_pt05_melbourne:       -999.0, $ ; Melbourne & Salzer (2002) procedure for pt05
;;;     zstrong_12oh_pt05_melbourne_err:   -999.0, $
;;;     zstrong_12oh_pt05_niiha:           -999.0, $ ; use [N II]/Ha to select the appropriate pt05 calibration, using the average 
;;;     zstrong_12oh_pt05_niiha_err:       -999.0, $ ;    of the two calibrations between 7.95 and 8.2
;;;     zstrong_ew_12oh_m91_o32:           -999.0, $ ; use ew(o32) to choose the branch
;;;     zstrong_ew_12oh_m91_o32_err:       -999.0, $
;;;      zstrong_12oh_oiiinii_niiha:        -999.0, $ ; use [N II]/Ha where ([O III]/Hb)/([N II]/Ha) is not measured
;;;      zstrong_12oh_oiiinii_niiha_err:    -999.0, $

;;; EW(N2)    
;;       
;;       if (tag_exist(line[0],'NII_6584_EW') and tag_exist(line[0],'H_ALPHA_EW')) then begin
;;          lineratio, line, 'NII_6584_EW', 'H_ALPHA_EW', '', '', $ ; NOTE!
;;            x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;;          if (nindex ne 0L) then begin
;;             abund[index].zstrong_ew_niiha = x
;;             abund[index].zstrong_ew_niiha_err = xerr
;;          endif
;;       endif

;;; ---------------------------------------------------------------------------    
;;; Zaritsky et al. (1994) - based on EWs
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_ew_r23 gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       x = abund[index].zstrong_ew_r23
;;       xerr = abund[index].zstrong_ew_r23_err
;;
;;       c = [9.265D,-0.33,-0.202,-0.207,-0.333]
;;       
;;       abund[index].zstrong_ew_12oh_zkh94 = c[0] + c[1]*x + c[2]*x^2 +c[3]*x^3 + c[4]*x^4
;;       abund[index].zstrong_ew_12oh_zkh94_err = sqrt( (c[1]*xerr)^2 + (2*c[2]*x*xerr)^2 + $
;;         (3*c[3]*x^2*xerr)^2 + (4*c[4]*x^3*xerr)^2 )
;;;      abund[index].zstrong_12oh_zkh94_err = 0.2
;;
;;       x = 0 & xerr = 0 ; save memory
;;       
;;; force the ZKH abundances to be greater than 8.4
;;
;;;      toosmall = where(abund[index].zstrong_12oh_zkh94 lt 8.4,ntoosmall)
;;;      if (ntoosmall ne 0L) then begin
;;;         abund[index[toosmall]].zstrong_12oh_zkh94 = -999.0
;;;         abund[index[toosmall]].zstrong_12oh_zkh94_err = -999.0
;;;      endif
;;       
;;    endif

;;; ---------------------------------------------------------------------------    
;;; use ew(o32) to select the appropriate branch of the m91 calibration;
;;; in the turn-around region use my empirical o32 calibration (based on
;;; the line-fluxes); define the turn-around region as being between 8.2
;;; and 8.5 (based on the o32-12+log(O/H) calibration
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_ew_o32 gt -900.0) and (abund.zstrong_ew_12oh_m91_upper gt -900.0) and $
;;      (abund.zstrong_ew_12oh_m91_lower gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       lower = where(abund[index].zstrong_ew_o32 gt 0.6,nlower)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_ew_12oh_m91_o32 = abund[index[lower]].zstrong_ew_12oh_m91_lower
;;          abund[index[lower]].zstrong_ew_12oh_m91_o32_err = abund[index[lower]].zstrong_ew_12oh_m91_lower_err
;;       endif
;;
;;       upper = where(abund[index].zstrong_ew_o32 lt 0.1,nupper)
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_ew_12oh_m91_o32 = abund[index[upper]].zstrong_ew_12oh_m91_upper
;;          abund[index[upper]].zstrong_ew_12oh_m91_o32_err = abund[index[upper]].zstrong_ew_12oh_m91_upper_err
;;       endif
;;
;;       turn = where((abund[index].zstrong_ew_o32 ge 0.1) and (abund[index].zstrong_ew_o32 le 0.6),nturn)
;;       if (nturn ne 0L) then begin
;;       
;;          c = [8.59D,-0.622D]
;;          cerr = [0.005D,0.010D]
;;
;;          x = abund[index[turn]].zstrong_ew_o32
;;          xerr = abund[index[turn]].zstrong_ew_o32_err
;;          
;;          abund[index[turn]].zstrong_ew_12oh_m91_o32 = c[0] + c[1] * x
;;          abund[index[turn]].zstrong_ew_12oh_m91_o32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;
;;       endif
;;
;;    endif

; if we need to normalize to H-beta then if H-beta is not
; well-measured then set that line to unmeasured as well; this should
; not be a big deal because we also cut on reddening, which requires
; H-beta; do not probagate the error in the normalization because that
; is done below
    
;   if keyword_set(normalize) then begin
;
;      tline = line
;      
;      linename = strtrim(tline[0].linename,2)
;      nline = n_elements(linename)
;      tags = tag_names(tline[0])
;
;      hbtag = where(strmatch(tags,'H_BETA',/fold) eq 1B,nmatch)
;      
;      for j = 0L, nspec-1L do begin
;
;         hbflux = reform((tline[j].(hbtag))[0,*])
;         hbferr = reform((tline[j].(hbtag))[1,*])
;
;         for k = 0L, nline-1L do begin 
;
;            match = where(strmatch(tags,linename[k],/fold) eq 1B,nmatch)
;            flux = (tline[j].(match))[0]
;            ferr = (tline[j].(match))[1]
;            if (ferr gt 0.0) and (hbferr gt 0.0) then begin
;               flux = flux / hbflux
;               ferr = ferr / hbflux
;               tline[j].(match) = transpose([ [flux], [ferr] ])
;            endif
;            
;         endfor
;
;      endfor
;
;      iline = tline[good_ebv
;   
;   endif else iline = line[good_ebv]
;
;   iabund = abund[good_ebv]

;;       index = where((abund.zstrong_ew_p gt -900.0) and (abund.zstrong_ew_r23 gt -900.0),nindex)
;;       if (nindex ne 0L) then begin
;;
;;; compute some ratios       
;;
;;          r23 = 10^abund[index].zstrong_ew_r23
;;          r23_err = alog(10.0)*abund[index].zstrong_ew_r23_err*r23
;;          r23_monte = rebin(reform(r23,1,nindex),nmonte,nindex)+randomn(seed,nmonte,nindex)*rebin(reform(r23_err,1,nindex),nmonte,nindex)
;;
;;          P = abund[index].zstrong_ew_P
;;          P_err = abund[index].zstrong_ew_P_err
;;          P_monte = rebin(reform(P,1,nindex),nmonte,nindex)+randomn(seed,nmonte,nindex)*rebin(reform(P_err,1,nindex),nmonte,nindex)
;;          
;;; lower branch - only works below 8.00
;;
;;          abund[index].zstrong_ew_12oh_pt05_lower = (r23 + 106.4 + 106.8*P - 3.40*P^2) / (17.72 + 6.60*P + 6.95*P^2 - 0.302*r23)
;;
;;          pt05_monte = (r23_monte + 106.4 + 106.8*P_monte - 3.40*P_monte^2) / $
;;            (17.72 + 6.60*P_monte + 6.95*P_monte^2 - 0.302*r23_monte)
;;          for ii = 0L, nindex-1L do abund[index[ii]].zstrong_ew_12oh_pt05_lower_err = stddev(pt05_monte[*,ii])
;;
;;          pt05_monte = 0        ; save memory
;;
;;; upper branch - only works above 8.25
;;
;;          abund[index].zstrong_ew_12oh_pt05_upper = (r23 + 726.1 + 842.2*P + 337.5*P^2) / (85.96 + 82.76*P + 43.98*P^2 + 1.793*r23)
;;
;;          pt05_monte = (r23_monte + 726.1 + 842.2*P_monte + 337.5*P_monte^2) / $
;;            (85.96 + 82.76*P_monte + 43.98*P_monte^2 + 1.793*r23_monte)
;;          for ii = 0L, nindex-1L do abund[index[ii]].zstrong_ew_12oh_pt05_upper_err = stddev(pt05_monte[*,ii])
;;
;;          pt05_monte = 0 & P_monte = 0 & r23_monte = 0 & r23 = 0 ; save memory
;;          r23_err = 0 & P = 0 & P_err = 0
;;          
;;       endif
;;; r3 (Pilyugin 2001)
;;    
;;    if (tag_exist(line[0],'OIII_4959') and tag_exist(line[0],'OIII_5007') and $
;;      tag_exist(line[0],'H_BETA')) then begin
;;
;;       lineratio, line, ['OIII_4959','OIII_5007'], 'H_BETA', '', '', $
;;         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;;
;;       if (nindex ne 0L) then begin
;;          abund[index].zstrong_r3 = x
;;          abund[index].zstrong_r3_err = xerr
;;       endif
;;
;;    endif
;;
;;; [O III]/H-beta
;;
;;    if (tag_exist(line[0],'H_BETA') and tag_exist(line[0],'OIII_5007')) then begin
;;
;;       lineratio, line, 'OIII_5007', 'H_BETA', '', '', $
;;         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
;;
;;       if (nindex ne 0L) then begin
;;          abund[index].zstrong_oiiihbeta = x
;;          abund[index].zstrong_oiiihbeta_err = xerr
;;       endif
;;
;;;;; ---------------------------------------------------------------------------    
;;;;; compute an empirical T(O++) based on [N II]/Ha then predict T(O+) 
;;;;; from the photionization models
;;;;; ---------------------------------------------------------------------------    
;;;;
;;;;    index = where(abund.zstrong_niiha gt -900.0,nindex)
;;;;    if (nindex ne 0L) then begin
;;;;
;;;;       c = [2512D,-7578D]
;;;;       cerr = [359D,239D]
;;;;
;;;;       x = abund[index].zstrong_niiha
;;;;       xerr = abund[index].zstrong_niiha_err
;;;;
;;;;       inrange = where((x gt -2.5) and (x lt -0.3),ninrange)
;;;;       if (ninrange ne 0L) then begin
;;;;          
;;;;          abund[index[inrange]].zstrong_T_oiii = c[0] + c[1] * x[inrange]
;;;;          abund[index[inrange]].zstrong_T_oiii_err = sqrt( cerr[0]^2 + (c[1] * xerr[inrange])^2.0 + $
;;;;            (cerr[1] * x[inrange])^2.0 )
;;;;
;;;;          t_oii = im_predict_toii(abund[index[inrange]].zstrong_T_oiii,toiii_lower=abund[index[inrange]].zstrong_T_oiii_err,$
;;;;            toiii_upper=abund[index[inrange]].zstrong_T_oiii_err)
;;;;
;;;;          abund[index[inrange]].zstrong_T_oii     = t_oii.ZT_T_oii
;;;;          abund[index[inrange]].zstrong_T_oii_err = t_oii.ZT_T_oii_err
;;;;          
;;;;          oh12Te_input = struct_addtags(line[index[inrange]],im_struct_trimtags(abund[index[inrange]],$
;;;;            select=['zstrong_DENSITY','zstrong_T_OIII','zstrong_T_OIII_ERR','zstrong_T_OII','zstrong_T_OII_ERR'],$
;;;;            newtags=['ZT_DENSITY','ZT_T_OIII','ZT_T_OIII_ERR','ZT_T_OII','ZT_T_OII_ERR']))
;;;;
;;;;          splog, 'Computing empirical electron temperature abundances.'
;;;;          oh12Te = im_compute_oh12te(oh12Te_input)
;;;;
;;;;          abund[index[inrange]].zstrong_log12oh     = oh12Te.ZT_log12oh
;;;;          abund[index[inrange]].zstrong_log12oh_err = oh12Te.ZT_log12oh_err
;;;;          abund[index[inrange]].zstrong_log_no      = oh12Te.ZT_log_no
;;;;          abund[index[inrange]].zstrong_log_no_err  = oh12Te.ZT_log_no_err
;;;;
;;;;       endif
;;;;          
;;;;    endif
;;;;
;;;; ---------------------------------------------------------------------------    
;;;; 12+log (O/H) based on my [N II]/Ha-Te calibration
;;;; ---------------------------------------------------------------------------    
;;;
;;;    index = where((abund.zstrong_niiha gt -900),nindex)
;;;;   index = where((abund.zstrong_niiha ge -2.5) and (abund.zstrong_niiha le -0.3),nindex)
;;;    if (nindex ne 0L) then begin
;;;
;;;       x = abund[index].zstrong_niiha
;;;       xerr = abund[index].zstrong_niiha_err
;;;
;;;       c = [9.153,0.782]
;;;       cerr = [0.0485,0.0307]
;;;
;;;;      abund[index].zstrong_12oh_niiha = c[0] + c[1] * x + c[2] * x^2
;;;;      abund[index].zstrong_12oh_niiha_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + $
;;;;        (cerr[1]*x)^2.0 + (2*c[2]*x*xerr)^2.0 + (cerr[2]*x^2)^2.0 )
;;;
;;;       abund[index].zstrong_12oh_niiha = c[0] + c[1] * x
;;;       abund[index].zstrong_12oh_niiha_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;;
;;;    endif
;;    
;;;;; ---------------------------------------------------------------------------    
;;;;; Moustakas et al. [N II]/Ha calibration; don't change the limits!
;;;;; ---------------------------------------------------------------------------    
;;;;
;;;;    index = where((abund.zstrong_niiha ge -2.5) and (abund.zstrong_niiha lt -0.3),nindex)
;;;;    if (nindex ne 0L) then begin
;;;;
;;;;       c = [9.193,0.7940]
;;;;       cerr = [0.04453,0.02749]
;;;;;      c = [9.2835823,0.85966437]
;;;;;      cerr = [0.058201499,0.036978448]
;;;;
;;;;       x = abund[index].zstrong_niiha
;;;;       xerr = abund[index].zstrong_niiha_err
;;;;
;;;;       abund[index].zstrong_12oh_niiha_moustakas = c[0] + c[1]*x
;;;;       abund[index].zstrong_12oh_niiha_moustakas_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;;;
;;;;    endif
;;    
;;;;; ---------------------------------------------------------------------------    
;;;;; Melbourne & Salzer (2002) [N II]/Ha calibration
;;;;; ---------------------------------------------------------------------------    
;;;;
;;;;    index = where(abund.zstrong_niiha gt -900.0,nindex)
;;;;    if (nindex ne 0L) then begin
;;;;
;;;;       c = [9.26D,1.23D,0.204D]
;;;;       cerr = [0.0D,0.0D,0.0D]
;;;;;      cerr = [0.156D,0.0D,0.0D]
;;;;
;;;;       x = abund[index].zstrong_niiha
;;;;       xerr = abund[index].zstrong_niiha_err
;;;;
;;;;       abund[index].zstrong_12oh_niiha_melbourne = c[0] + c[1]*x + c[2]*x^2
;;;;       abund[index].zstrong_12oh_niiha_melbourne_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 + $
;;;;         (2*c[2]*x*xerr)^2 + (cerr[2]*x^2)^2 )
;;;;
;;;;    endif
;;    
;;;; ---------------------------------------------------------------------------    
;;;; Denicolo et al. (2002) [N II]/Ha calibration
;;;; ---------------------------------------------------------------------------    
;;;
;;;    index = where(abund.zstrong_niiha gt -900.0,nindex)
;;;    if (nindex ne 0L) then begin
;;;
;;;       c = [9.12D,0.73D]
;;;       cerr = [0.05D,0.10D]
;;;       
;;;       x = abund[index].zstrong_niiha
;;;       xerr = abund[index].zstrong_niiha_err
;;;       
;;;       abund[index].zstrong_12oh_niiha_denicolo = c[0] + c[1] * x
;;;       abund[index].zstrong_12oh_niiha_denicolo_err = sqrt( cerr[0]^2 + (c[1] * xerr)^2.0 + (cerr[1] * x)^2.0 )
;;;
;;;       x = 0 & xerr = 0 ; save memory
;;;
;;;    endif
    
;; ---------------------------------------------------------------------------    
;; 12+log (O/H) based on my o32-Te calibration
;; ---------------------------------------------------------------------------    
;
;    index = where((abund.zstrong_o32 ge -0.9) and (abund.zstrong_o32 le 1.2),nindex)
;    if (nindex ne 0L) then begin
;
;       c = [8.2835434,-0.47751407]
;       cerr = [0.024664146,0.042504591]
;
;       x = abund[index].zstrong_o32
;       xerr = abund[index].zstrong_o32_err
;
;       abund[index].zstrong_12oh_o32 = c[0] + c[1] * x
;       abund[index].zstrong_12oh_o32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;
;    endif
    
;; ---------------------------------------------------------------------------    
;; 12+log (O/H) based on my P-Te calibration
;; ---------------------------------------------------------------------------    
;
;    index = where((abund.zstrong_P ge 0.1) and (abund.zstrong_P le 1.0),nindex)
;;   index = where((abund.zstrong_P ge -900.0),nindex)
;    if (nindex ne 0L) then begin
;
;       c = [8.8139881,-1.0486419]
;       cerr = [0.067280689,0.095832449]
;
;       x = abund[index].zstrong_P
;       xerr = abund[index].zstrong_P_err
;
;       abund[index].zstrong_12oh_P = c[0] + c[1] * x
;       abund[index].zstrong_12oh_P_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;
;    endif
    
;;; ---------------------------------------------------------------------------    
;;; 12+log (O/H) based on my ([O III]/Hb)/([N II]/Ha)-Te calibration - WRONG!!!
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_oiiinii ge -0.7) and (abund.zstrong_oiiinii le 1.9),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       c = [8.5638200,-0.21535851]
;;       cerr = [0.046192925,0.032765583]
;;
;;       x = abund[index].zstrong_oiiinii
;;       xerr = abund[index].zstrong_oiiinii_err
;;
;;       abund[index].zstrong_12oh_oiiinii = c[0] + c[1] * x
;;       abund[index].zstrong_12oh_oiiinii_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;
;;    endif
    
;;; ---------------------------------------------------------------------------    
;;; Moustakas et al. ([O III]/Hb)/([N II]/Ha) calibration
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_oiiinii ge -1.0) and (abund.zstrong_oiiinii le 1.9),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       c = [8.566,-0.2244]
;;       cerr = [0.04537,0.031538]
;;;      c = [8.5729971,-0.24480717]
;;;      cerr = [0.055266465,0.038004130]
;;
;;       x = abund[index].zstrong_oiiinii
;;       xerr = abund[index].zstrong_oiiinii_err
;;
;;       abund[index].zstrong_12oh_oiiinii_moustakas = c[0] + c[1] * x
;;       abund[index].zstrong_12oh_oiiinii_moustakas_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;
;;    endif
    
;;; ---------------------------------------------------------------------------    
;;; Charlot & Longhetti (2001, Case F)    
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_oiiihbeta gt -900.0) and (abund.zstrong_oiioiii gt -900.0),nindex)
;;
;;    if (nindex ne 0L) then begin
;;
;;       logoiioiii = abund[index].zstrong_oiioiii
;;       logoiioiii_err = abund[index].zstrong_oiioiii_err       
;;
;;       oiioiii = 10^logoiioiii
;;       oiioiii_err = alog(10.0)*logoiioiii_err*oiioiii
;;       
;;       x2 = oiioiii/1.5
;;       x2_err = oiioiii_err/1.5
;;
;;       logoiiihbeta = abund[index].zstrong_oiiihbeta       
;;       logoiiihbeta_err = abund[index].zstrong_oiiihbeta_err
;;
;;       oiiihbeta = 10^logoiiihbeta
;;       oiiihbeta_err = alog(10.0)*logoiiihbeta_err*oiiihbeta
;;       
;;       x3 = oiiihbeta/2.0
;;       x3_err = oiiihbeta_err/2.0
;;       
;;       lo = where(oiioiii lt 0.8,nlo,comp=hi,ncomp=nhi)
;;       if (nlo ne 0L) then begin
;;
;;          term = 3.78D-4 * x2[lo]^(0.17) * x3[lo]^(-0.44)
;;          term_err = 3.78D-4 * sqrt( (0.17*x2[lo]^(0.17-1)*x3[lo]^(-0.44)*x2_err[lo])^2 + $
;;            (-0.44*x3[lo]^(-0.44-1)*x2[lo]^0.17*x3_err[lo])^2 )
;;
;;          abund[index[lo]].zstrong_12oh_cl01F = 12 + alog10(term)
;;          abund[index[lo]].zstrong_12oh_cl01F_err = term_err/term/alog(10.0)
;;;         abund[index[lo]].zstrong_12oh_cl01F_err = 0.51
;;
;;       endif
;;       
;;       if (nhi ne 0L) then begin
;;
;;          term = 3.96D-4 * x3[hi]^(-0.46)
;;          term_err = sqrt( (3.96D-4*0.46*x3[hi]^(-0.46-1)*x3_err[hi])^2 )
;;          
;;          abund[index[hi]].zstrong_12oh_cl01F = 12 + alog10(term)
;;          abund[index[hi]].zstrong_12oh_cl01F_err = term_err/term/alog(10.0)
;;;         abund[index[hi]].zstrong_12oh_cl01F_err = 0.31
;;
;;       endif
;;
;;    endif
    
;;; ---------------------------------------------------------------------------    
;;; use the procedure in Contini et al. (2002), Section 4.1, to
;;; determine if an object is in the upper or lower branch of the m91
;;; relation
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_niiha gt -900.0) and (abund.zstrong_niioii gt -900.0) and $
;;      (abund.zstrong_12oh_m91_upper gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;; [N II]/Ha and [N II]/[O II] both indicate the lower branch
;;       
;;       lower = where((abund[index].zstrong_niiha lt -1.0) and (abund[index].zstrong_niioii lt -1.05),nlower)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_12oh_m91_contini = abund[index[lower]].zstrong_12oh_m91_lower
;;          abund[index[lower]].zstrong_12oh_m91_contini_err = abund[index[lower]].zstrong_12oh_m91_lower_err
;;       endif
;;
;;; [N II]/Ha and [N II]/[O II] both indicate the upper branch
;;
;;       upper = where((abund[index].zstrong_niiha gt -1.0) and (abund[index].zstrong_niioii gt -0.8),nupper)
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_12oh_m91_contini = abund[index[upper]].zstrong_12oh_m91_upper
;;          abund[index[upper]].zstrong_12oh_m91_contini_err = abund[index[upper]].zstrong_12oh_m91_upper_err
;;       endif
;;
;;; [N II]/[O II] is in an ambigious region --> use [N II]/Ha
;;
;;       ambig = where((abund[index].zstrong_niioii gt -1.05) and (abund[index].zstrong_niioii lt -0.8),nambig)
;;
;;       if (nambig ne 0L) then begin
;;
;;          lower = where((abund[index[ambig]].zstrong_niiha lt -1.0),nlower,comp=upper,ncomp=nupper)
;;
;;          if (nlower ne 0L) then begin
;;             abund[index[ambig[lower]]].zstrong_12oh_m91_contini = abund[index[ambig[lower]]].zstrong_12oh_m91_lower
;;             abund[index[ambig[lower]]].zstrong_12oh_m91_contini_err = abund[index[ambig[lower]]].zstrong_12oh_m91_lower_err
;;          endif
;;
;;          if (nupper ne 0L) then begin
;;             abund[index[ambig[upper]]].zstrong_12oh_m91_contini = abund[index[ambig[upper]]].zstrong_12oh_m91_upper
;;             abund[index[ambig[upper]]].zstrong_12oh_m91_contini_err = abund[index[ambig[upper]]].zstrong_12oh_m91_upper_err
;;          endif
;;
;;       endif
;;
;;; [N II]/Ha indicates upper branch but [N II]/[O II] indicates lower
;;; branch --> take an average of the two metallicities
;;       
;;       veryambig = where((abund[index].zstrong_niiha gt -1.0) and (abund[index].zstrong_niioii lt -1.05),nambig)
;;
;;       if (nambig ne 0L) then begin
;;
;;          Z1 = abund[index[veryambig]].zstrong_12oh_m91_lower
;;          Z1err = abund[index[veryambig]].zstrong_12oh_m91_lower_err
;;
;;          Z2 = abund[index[veryambig]].zstrong_12oh_m91_upper
;;          Z2err = abund[index[veryambig]].zstrong_12oh_m91_upper_err
;;          
;;          abund[index[veryambig]].zstrong_12oh_m91_contini = (Z1 + Z2) / 2.0
;;          abund[index[veryambig]].zstrong_12oh_m91_contini_err = sqrt(Z1err^2 + Z2err^2)
;;
;;       endif
;;          
;;; [N II]/Ha indicates lower branch but [N II]/[O II] indicates upper
;;; branch --> take an average of the two metallicities (should rarely
;;; occur) 
;;       
;;       veryambig = where((abund[index].zstrong_niiha lt -1.0) and (abund[index].zstrong_niioii gt -0.8),nambig)
;;
;;       if (nambig ne 0L) then begin
;;
;;          Z1 = abund[index[veryambig]].zstrong_12oh_m91_lower
;;          Z1err = abund[index[veryambig]].zstrong_12oh_m91_lower_err
;;
;;          Z2 = abund[index[veryambig]].zstrong_12oh_m91_upper
;;          Z2err = abund[index[veryambig]].zstrong_12oh_m91_upper_err
;;          
;;          abund[index[veryambig]].zstrong_12oh_m91_contini = (Z1 + Z2) / 2.0
;;          abund[index[veryambig]].zstrong_12oh_m91_contini_err = sqrt(Z1err^2 + Z2err^2)
;;
;;       endif
;;
;;    endif

;;; ---------------------------------------------------------------------------    
;;; use zkh94 or cl01A to select the appropriate branch according to
;;; Kewley & Dopita (2002)
;;; ---------------------------------------------------------------------------    
;;
;;;   index = where((abund.zstrong_12oh_cl01A gt -900.0) and (abund.zstrong_12oh_m91_upper gt -900.0),nindex)
;;    index = where((abund.zstrong_12oh_zkh94 gt -900.0) and (abund.zstrong_12oh_m91_upper gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;;      lower = where(abund[index].zstrong_12oh_cl01A le 8.4,nlower,comp=upper,ncomp=nupper)
;;       lower = where(abund[index].zstrong_12oh_zkh94 le 8.4,nlower,comp=upper,ncomp=nupper)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_12oh_m91_kd02 = abund[index[lower]].zstrong_12oh_m91_lower
;;          abund[index[lower]].zstrong_12oh_m91_kd02_err = abund[index[lower]].zstrong_12oh_m91_lower_err
;;       endif
;;
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_12oh_m91_kd02 = abund[index[upper]].zstrong_12oh_m91_upper
;;          abund[index[upper]].zstrong_12oh_m91_kd02_err = abund[index[upper]].zstrong_12oh_m91_upper_err
;;       endif
;;       
;;    endif
    
;;; ---------------------------------------------------------------------------    
;;; use o32 to select the appropriate branch of the m91 calibration; in
;;; the turn-around region use my empirical o32 calibration, but add the 
;;; 0.20+/-0.05 systematic offset between the two abundance scales
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_o32 gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       lower = where((abund[index].zstrong_o32 gt 0.55) and (abund[index].zstrong_12oh_m91_lower gt -900.0),nlower)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_12oh_m91_o32     = abund[index[lower]].zstrong_12oh_m91_lower
;;          abund[index[lower]].zstrong_12oh_m91_o32_err = abund[index[lower]].zstrong_12oh_m91_lower_err
;;       endif
;;
;;       upper = where((abund[index].zstrong_o32 lt 0.11) and (abund[index].zstrong_12oh_m91_upper gt -900.0),nupper)
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_12oh_m91_o32     = abund[index[upper]].zstrong_12oh_m91_upper
;;          abund[index[upper]].zstrong_12oh_m91_o32_err = abund[index[upper]].zstrong_12oh_m91_upper_err
;;       endif
;;
;;       turn = where((abund[index].zstrong_o32 ge 0.11) and (abund[index].zstrong_o32 le 0.55),nturn)
;;       if (nturn ne 0L) then begin
;;
;;; bisector fit to Atlas+NFGS+HII Regions; the scatter perpendicular
;;; to the bisector is 0.08 dex = systematic error
;;
;;          offset = 0.20
;;          offset_err = 0.05
;;       
;;          sys_err = 0.08
;;          
;;;         c = [8.329D,-0.512D]
;;          c = [8.328D,-0.473D]
;;          cerr = [0.005D,0.013D]
;;
;;          x = abund[index[turn]].zstrong_o32
;;          xerr = abund[index[turn]].zstrong_o32_err
;;          
;;          abund[index[turn]].zstrong_12oh_m91_o32 = c[0] + c[1] * x + offset
;;          abund[index[turn]].zstrong_12oh_m91_o32_err = sqrt( cerr[0]^2 + $
;;            (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 + offset_err^2 + sys_err^2)
;;
;;       endif
;;
;;    endif

;;; ---------------------------------------------------------------------------    
;;; use [O III]/Hb to select the appropriate branch of the kk04 calibration
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_oiiihbeta gt -900.0) and (abund.zstrong_12oh_kk04_r23_upper gt -900.0) and $
;;      (abund.zstrong_12oh_kk04_r23_lower gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       lower = where(abund[index].zstrong_oiiihbeta gt 0.25,nlower,comp=upper,ncomp=nupper)
;;       if (nlower ne 0L) then abund[index[lower]].zstrong_12oh_kk04_o32 = abund[index[lower]].zstrong_12oh_kk04_r23_lower
;;       if (nupper ne 0L) then abund[index[upper]].zstrong_12oh_kk04_o32 = abund[index[upper]].zstrong_12oh_kk04_r23_upper
;;       
;;    endif
    
; ---------------------------------------------------------------------------    
; Pilyugin (2000, 2001) - old calibration
; ---------------------------------------------------------------------------    

;   index = where((abund.zstrong_r3 gt -900.0) and (abund.zstrong_r23 gt -900.0),nindex)
;   if (nindex ne 0L) then begin
;
; compute some ratios       
;
;      logr23 = abund[index].zstrong_r23
;      logr23_err = abund[index].zstrong_r23_err
;      
;      r23 = 10^logr23
;      r23_err = alog(10.0)*logr23_err*r23
;
;      logr3 = abund[index].zstrong_r3
;      logr3_err = abund[index].zstrong_r3_err
;      
;      r3 = 10^logr3
;      r3_err = alog(10.0)*logr3_err*r3
;      
;      P = r3/r23
;      P_err = sqrt( (r3_err/r23)^2 + (r23_err*r3/r23^2)^2 )
;
; lower branch
;
;      abund[index].zstrong_12oh_p01_lower = 6.35 + 3.19*logr23 - 1.74*logr3
;      abund[index].zstrong_12oh_p01_lower_err = sqrt( (3.19*logr23_err)^2 + (1.74*logr3_err)^2 )
;
; upper branch
;
;      top = r23 + 54.2 + 59.45*P + 7.31*P^2
;      top_err = sqrt( r23_err^2 + (59.45*P_err)^2 + (2.0*7.31*P*P_err)^2 )
;
;      but = 6.07 + 6.71*P + 0.37*P^2 + 0.243*r23
;      but_err = sqrt( (6.71*P_err)^2  + (2.0*0.37*P*P_err)^2 + (0.243*r23_err)^2 )
;
;      p01 = top/but
;      p01_err = sqrt( (top_err/but)^2 + (but_err*top/but^2)^2 )
;      
;      abund[index].zstrong_12oh_p01_upper = p01
;      abund[index].zstrong_12oh_p01_upper_err = p01_err

; compute the average of the two calibrations, to be used below

;      abund[index].zstrong_12oh_p01_average = (abund[index].zstrong_12oh_p01_lower + abund[index].zstrong_12oh_p01_upper) / 2.0
;      abund[index].zstrong_12oh_p01_average_err = sqrt((abund[index].zstrong_12oh_p01_lower_err^2 + abund[index].zstrong_12oh_p01_upper_err^2) / 2.0 )

; now, constrain the lower-branch abundances to be less than 7.95

;      toobig = where(abund[index].zstrong_12oh_p01_lower gt 7.95,ntoobig)
;      if (ntoobig ne 0L) then begin
;         abund[index[toobig]].zstrong_12oh_p01_lower = -999.0
;         abund[index[toobig]].zstrong_12oh_p01_lower_err = -999.0
;      endif
       
; and constrain the upper-branch abundances to be greater than 8.2

;      toosmall = where(abund[index].zstrong_12oh_p01_upper lt 8.2,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         abund[index[toosmall]].zstrong_12oh_p01_upper = -999.0
;         abund[index[toosmall]].zstrong_12oh_p01_upper_err = -999.0
;      endif
       
;   endif
       
; ---------------------------------------------------------------------------
; use the procedure outlined in Melbourne & Salzer (2002) to figure
; out which Pilyugin branch we're on but use the new upper-branch
; calibration by Pilyugin (2002) instead of the Edmunds & Pagel
; calibration used by Melbourne & Salzer
; ---------------------------------------------------------------------------

;;    index = where((abund.zstrong_niiha gt -900.0) and (abund.zstrong_12oh_pt05_upper gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;; [N II/Ha < -1.3 --> lower branch
;;       
;;       lower = where((abund[index].zstrong_niiha lt -1.3),nlower)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_12oh_pt05_melbourne = abund[index[lower]].zstrong_12oh_pt05_lower
;;          abund[index[lower]].zstrong_12oh_pt05_melbourne_err = abund[index[lower]].zstrong_12oh_pt05_lower_err
;;       endif
;;       
;;; [N II/Ha > -1.0 --> upper branch
;;       
;;       upper = where((abund[index].zstrong_niiha gt -1.0),nupper)
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_12oh_pt05_melbourne = abund[index[upper]].zstrong_12oh_pt05_upper
;;          abund[index[upper]].zstrong_12oh_pt05_melbourne_err = abund[index[upper]].zstrong_12oh_pt05_upper_err
;;       endif
;;       
;;    endif

;;; ---------------------------------------------------------------------------    
;;; use o32 to select the appropriate branch of the pt05 calibration; in
;;; the turn-around region use my empirical o32 calibration
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((abund.zstrong_o32 gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       lower = where((abund[index].zstrong_o32 gt 0.749) and (abund[index].zstrong_12oh_pt05_lower gt -900.0),nlower)
;;       if (nlower ne 0L) then begin
;;          abund[index[lower]].zstrong_12oh_pt05_o32     = abund[index[lower]].zstrong_12oh_pt05_lower
;;          abund[index[lower]].zstrong_12oh_pt05_o32_err = abund[index[lower]].zstrong_12oh_pt05_lower_err
;;       endif
;;
;;       upper = where((abund[index].zstrong_o32 lt 0.308) and (abund[index].zstrong_12oh_pt05_upper gt -900.0),nupper)
;;       if (nupper ne 0L) then begin
;;          abund[index[upper]].zstrong_12oh_pt05_o32     = abund[index[upper]].zstrong_12oh_pt05_upper
;;          abund[index[upper]].zstrong_12oh_pt05_o32_err = abund[index[upper]].zstrong_12oh_pt05_upper_err
;;       endif
;;
;;       turn = where((abund[index].zstrong_o32 ge 0.308) and (abund[index].zstrong_o32 le 0.749),nturn)
;;       if (nturn ne 0L) then begin
;;
;;; bisector fit to Atlas+NFGS+HII Regions; the scatter perpendicular
;;; to the bisector is 0.08 dex = systematic error
;;
;;          sys_err = 0.08
;;          
;;;         c = [8.329D,-0.512D]
;;          c = [8.328D,-0.473D]
;;          cerr = [0.005D,0.013D]
;;
;;          x = abund[index[turn]].zstrong_o32
;;          xerr = abund[index[turn]].zstrong_o32_err
;;          
;;          abund[index[turn]].zstrong_12oh_pt05_o32 = c[0] + c[1] * x
;;          abund[index[turn]].zstrong_12oh_pt05_o32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 + sys_err^2)
;;
;;       endif
;;
;;    endif

;;; ---------------------------------------------------------------------------
;;; use [N II]/Ha where ([O III]/Hb)/([N II]/Ha) is not measured and
;;; where [N II]/Ha implies a metallicity less than 8.5
;;; ---------------------------------------------------------------------------
;;
;;    index = where((abund.zstrong_12oh_oiiinii_pettini gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       abund[index].zstrong_12oh_oiiinii_niiha = abund[index].zstrong_12oh_oiiinii_pettini
;;       abund[index].zstrong_12oh_oiiinii_niiha_err = abund[index].zstrong_12oh_oiiinii_pettini_err
;;
;;    endif
;;
;;    index = where((abund.zstrong_12oh_oiiinii_pettini lt -900.0) and (abund.zstrong_12oh_niiha_pettini gt -900.0) and $
;;      (abund.zstrong_12oh_niiha_pettini lt 8.5),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       abund[index].zstrong_12oh_oiiinii_niiha = abund[index].zstrong_12oh_niiha_pettini
;;       abund[index].zstrong_12oh_oiiinii_niiha_err = abund[index].zstrong_12oh_niiha_pettini_err
;;
;;    endif
    
;; ---------------------------------------------------------------------------
;; use [N II]/Ha [Pettini] to select the appropriate pt05 calibration; we could
;; use the average of the two calibrations between 12+log(O/H)=7.95 and
;; 8.2 (and the commented code below allows you to do that), but this
;; does not work well; so ignore the region between 7.95 and 8.2
;; ---------------------------------------------------------------------------
;
;    index = where((abund.zstrong_12oh_niiha_pettini gt -900.0) and ((abund.zstrong_12oh_pt05_upper gt -900.0) or $
;      (abund.zstrong_12oh_pt05_lower gt -900.0)),nindex)
;    if (nindex ne 0L) then begin
;
;       lower = where((abund[index].zstrong_12oh_niiha_pettini lt 7.95) and (abund[index].zstrong_12oh_pt05_lower gt -900),nlower)
;       if (nlower ne 0L) then begin
;          abund[index[lower]].zstrong_12oh_pt05_niiha = abund[index[lower]].zstrong_12oh_pt05_lower
;          abund[index[lower]].zstrong_12oh_pt05_niiha_err = abund[index[lower]].zstrong_12oh_pt05_lower_err
;       endif
;       
;       upper = where((abund[index].zstrong_12oh_niiha_pettini gt 8.2) and (abund[index].zstrong_12oh_pt05_upper gt -900),nupper)
;       if (nupper ne 0L) then begin
;          abund[index[upper]].zstrong_12oh_pt05_niiha = abund[index[upper]].zstrong_12oh_pt05_upper
;          abund[index[upper]].zstrong_12oh_pt05_niiha_err = abund[index[upper]].zstrong_12oh_pt05_upper_err
;       endif
;       
;;      turn = where((abund[index].zstrong_12oh_niiha_pettini ge 7.95) and (abund[index].zstrong_12oh_niiha_pettini le 8.2) and $
;;        (abund[index].zstrong_12oh_pt05_average gt -900),nturn)
;;      if(nturn ne 0L) then begin
;;         abund[index[turn]].zstrong_12oh_pt05_niiha = abund[index[turn]].zstrong_12oh_pt05_average
;;         abund[index[turn]].zstrong_12oh_pt05_niiha_err = abund[index[turn]].zstrong_12oh_pt05_average_err
;;      endif
;
;    endif

;;;; ---------------------------------------------------------------------------
;;;; use pt05 on the upper and lower branches but the Melbourne & Salzer
;;;; (2002) calibration of [N II]/Ha in the turn-around region (see van
;;;; Zee et al. 1998)
;;;; ---------------------------------------------------------------------------
;;;
;;;    index = where((abund.zstrong_12oh_niiha gt -900.0) or (abund.zstrong_12oh_P gt -900.0),nindex)
;;;;   index = where((abund.zstrong_12oh_niiha gt -900.0) or (abund.zstrong_12oh_P gt -900.0) and $
;;;;     (abund.zstrong_12oh_pt05_upper gt -900.0),nindex)
;;;    if (nindex ne 0L) then begin
;;;
;;;;;; if both upper & lower branch pt05 calibrations give the same branch
;;;;;; then use them!
;;;;;
;;;;;       upper = where((abund[index].zstrong_12oh_pt05_lower ge 8.2) and (abund[index].zstrong_12oh_pt05_upper ge 8.2),nupper)
;;;;;       if (nupper ne 0L) then begin
;;;;;          abund[index[upper]].zstrong_12oh_pt05 = abund[index[upper]].zstrong_12oh_pt05_upper
;;;;;          abund[index[upper]].zstrong_12oh_pt05_err = abund[index[upper]].zstrong_12oh_pt05_upper_err
;;;;;       endif
;;;;;       
;;;;;       lower = where((abund[index].zstrong_12oh_pt05_lower le 7.95) and (abund[index].zstrong_12oh_pt05_upper le 7.95),nlower)
;;;;;       if (nlower ne 0L) then begin
;;;;;          abund[index[lower]].zstrong_12oh_pt05 = abund[index[lower]].zstrong_12oh_pt05_lower
;;;;;          abund[index[lower]].zstrong_12oh_pt05_err = abund[index[lower]].zstrong_12oh_pt05_lower_err
;;;;;       endif
;;;
;;;; if we have [N II]/Ha and P then use [N II]/Ha otherwise use P
;;;
;;;       n2go = where((abund[index].zstrong_12oh_niiha gt -900.0),nn2go)
;;;;      n2go = where((abund[index].zstrong_12oh_niiha gt -900.0) and (abund[index].zstrong_12oh_pt05 lt -900.0),nn2go)
;;;       if (nn2go ne 0L) then begin
;;;
;;;; upper branch          
;;;          
;;;          upper = where((abund[index[n2go]].zstrong_12oh_niiha ge 8.2) and $
;;;            (abund[index[n2go]].zstrong_12oh_pt05_upper gt -900.0),nupper)
;;;;         upper = where(abund[index[n2go]].zstrong_12oh_niiha ge 8.2,nupper)
;;;
;;;          if (nupper ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nupper,format='(I0)')+$
;;;               ' upper-branch objects using niiha-oh.'
;;;             abund[index[n2go[upper]]].zstrong_12oh_pt05 = abund[index[n2go[upper]]].zstrong_12oh_pt05_upper
;;;             abund[index[n2go[upper]]].zstrong_12oh_pt05_err = abund[index[n2go[upper]]].zstrong_12oh_pt05_upper_err
;;;          endif
;;;          
;;;; lower branch          
;;;
;;;          lower = where((abund[index[n2go]].zstrong_12oh_niiha le 7.95) and $
;;;            (abund[index[n2go]].zstrong_12oh_pt05_lower gt -900.0),nlower)
;;;;         lower = where(abund[index[n2go]].zstrong_12oh_niiha le 7.95,nlower)
;;;
;;;          if (nlower ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nlower,format='(I0)')+$
;;;               ' lower-branch objects using [N II]/Ha.'
;;;             abund[index[n2go[lower]]].zstrong_12oh_pt05 = abund[index[n2go[lower]]].zstrong_12oh_pt05_lower
;;;             abund[index[n2go[lower]]].zstrong_12oh_pt05_err = abund[index[n2go[lower]]].zstrong_12oh_pt05_lower_err
;;;          endif
;;;
;;;; ambigious --> use [N II]/Ha 
;;;
;;;          ambig = where((abund[index[n2go]].zstrong_12oh_niiha gt 7.95) and (abund[index[n2go]].zstrong_12oh_niiha lt 8.2),nambig)
;;;
;;;          if (nambig ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Applying [N II]/Ha on '+string(nambig,format='(I0)')+' ambigious objects.'
;;;             abund[index[n2go[ambig]]].zstrong_12oh_pt05 = abund[index[n2go[ambig]]].zstrong_12oh_niiha
;;;             abund[index[n2go[ambig]]].zstrong_12oh_pt05_err = abund[index[n2go[ambig]]].zstrong_12oh_niiha_err
;;;          endif
;;;
;;;;;; ambigious --> use the lower/upper average
;;;;;
;;;;;          ambig = where((abund[index[n2go]].zstrong_12oh_niiha gt 7.95) and (abund[index[n2go]].zstrong_12oh_niiha lt 8.2),nambig)
;;;;;          if (nambig ne 0L) then begin
;;;;;             abund[index[n2go[ambig]]].zstrong_12oh_pt05 = $
;;;;;               (abund[index[n2go[ambig]]].zstrong_12oh_pt05_lower + abund[index[n2go[ambig]]].zstrong_12oh_pt05_upper) / 2.0
;;;;;             abund[index[n2go[ambig]]].zstrong_12oh_pt05_err = sqrt( abund[index[n2go[ambig]]].zstrong_12oh_pt05_lower_err^2 + $
;;;;;               abund[index[n2go[ambig]]].zstrong_12oh_pt05_upper_err^2)
;;;;;          endif
;;;
;;;       endif
;;;
;;;       Pgo = where((abund[index].zstrong_12oh_P gt -900.0) and (abund[index].zstrong_12oh_pt05 lt -900.0),nPgo)
;;;       if (nPgo ne 0L) then begin
;;;
;;;          upper = where((abund[index[Pgo]].zstrong_12oh_P ge 8.2) and $
;;;            (abund[index[Pgo]].zstrong_12oh_pt05_upper gt -900.0),nupper)
;;;
;;;          if (nupper ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nupper,format='(I0)')+$
;;;               ' upper-branch objects using P-oh.'
;;;             abund[index[Pgo[upper]]].zstrong_12oh_pt05 = abund[index[Pgo[upper]]].zstrong_12oh_pt05_upper
;;;             abund[index[Pgo[upper]]].zstrong_12oh_pt05_err = abund[index[Pgo[upper]]].zstrong_12oh_pt05_upper_err
;;;          endif
;;;          
;;;          lower = where((abund[index[Pgo]].zstrong_12oh_P le 7.95) and $
;;;            (abund[index[Pgo]].zstrong_12oh_pt05_lower gt -900.0),nlower)
;;;
;;;          if (nlower ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nlower,format='(I0)')+$
;;;               ' lower-branch objects using P-oh.'
;;;             abund[index[Pgo[lower]]].zstrong_12oh_pt05 = abund[index[Pgo[lower]]].zstrong_12oh_pt05_lower
;;;             abund[index[Pgo[lower]]].zstrong_12oh_pt05_err = abund[index[Pgo[lower]]].zstrong_12oh_pt05_lower_err
;;;          endif
;;;
;;;; ambigious --> use P
;;;
;;;          ambig = where((abund[index[Pgo]].zstrong_12oh_P gt 7.95) and (abund[index[Pgo]].zstrong_12oh_P lt 8.2),nambig)
;;;          if (nambig ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Applying P-oh on '+string(nambig,format='(I0)')+' ambigious objects.'
;;;             abund[index[Pgo[ambig]]].zstrong_12oh_pt05 = abund[index[Pgo[ambig]]].zstrong_12oh_P
;;;             abund[index[Pgo[ambig]]].zstrong_12oh_pt05_err = abund[index[Pgo[ambig]]].zstrong_12oh_P_err
;;;          endif
;;;
;;;;;; ambigious --> use the lower/upper average
;;;;;
;;;;;          ambig = where((abund[index[Pgo]].zstrong_12oh_P gt 7.95) and (abund[index[Pgo]].zstrong_12oh_P lt 8.2),nambig)
;;;;;          if (nambig ne 0L) then begin
;;;;;             abund[index[Pgo[ambig]]].zstrong_12oh_pt05 = $
;;;;;               (abund[index[Pgo[ambig]]].zstrong_12oh_pt05_lower + abund[index[Pgo[ambig]]].zstrong_12oh_pt05_upper) / 2.0
;;;;;             abund[index[Pgo[ambig]]].zstrong_12oh_pt05_err = sqrt( abund[index[Pgo[ambig]]].zstrong_12oh_pt05_lower_err^2 + $
;;;;;               abund[index[Pgo[ambig]]].zstrong_12oh_pt05_upper_err^2)
;;;;;          endif
;;;
;;;       endif
;;;       
;;;;      upper = where(abund[index].zstrong_niioii gt -0.7,nupper)
;;;;      if (nupper ne 0L) then abund[index[upper]].zstrong_12oh_pt05 = abund[index[upper]].zstrong_12oh_pt05_upper
;;;;      
;;;;      lower = where(abund[index].zstrong_niioii lt -1.3,nlower)
;;;;      if (nlower ne 0L) then abund[index[lower]].zstrong_12oh_pt05 = abund[index[lower]].zstrong_12oh_pt05_lower
;;;;
;;;;      ambig = where(abund[index].zstrong_12oh_pt05 lt -900.0,nambig)
;;;;      if (nambig ne 0L) then begin
;;;;      
;;;;         upper_diff = abs(abund[index[ambig]].zstrong_12oh_niiha - abund[index[ambig]].zstrong_12oh_pt05_upper)
;;;;         lower_diff = abs(abund[index[ambig]].zstrong_12oh_niiha - abund[index[ambig]].zstrong_12oh_pt05_lower)
;;;;
;;;;         toobig = where((upper_diff ge 0.05) and (lower_diff ge 0.05),ntoobig,comp=okay,ncomp=nokay)
;;;;         if (ntoobig ne 0L) then abund[index[ambig[toobig]]].zstrong_12oh_pt05 = $
;;;;           abund[index[ambig[toobig]]].zstrong_12oh_niiha
;;;;
;;;;         if (nokay ne 0L) then begin
;;;;         
;;;;            upper = where((upper_diff[okay] lt lower_diff[okay]),nupper)
;;;;            if (nupper ne 0L) then abund[index[ambig[okay[upper]]]].zstrong_12oh_pt05 = $
;;;;              abund[index[ambig[okay[upper]]]].zstrong_12oh_pt05_upper
;;;;         
;;;;            lower = where((lower_diff[okay] lt upper_diff[okay]),nlower)
;;;;            if (nlower ne 0L) then abund[index[ambig[okay[lower]]]].zstrong_12oh_pt05 = $
;;;;              abund[index[ambig[okay[lower]]]].zstrong_12oh_pt05_lower
;;;;
;;;;         endif
;;;;
;;;;      endif
;;;
;;;    endif 

