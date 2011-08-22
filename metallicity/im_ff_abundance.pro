;+
; NAME:
;   IM_FF_ABUNDANCE()
;
; PURPOSE:
;   Compute abundances using Pilyugin's "ff" method.
;
; INPUTS:
;   line   - iSPEC-style input structure of de-reddened line fluxes 
;
; OPTIONAL INPUTS:
;   snrcut_abundance - require that individual emission lines have
;     a signal-to-noise ratio greater than SNRCUT_ABUNDANCE 
;   nmonte - number of Monte Carlo realizations for computing the
;     errors in particular abundances (default 500)
;
; KEYWORD PARAMETERS:
;   electrondensity - compute the electron density using the sulfur
;     line-ratio 
;   observed - use the observed rather than the predicted value of R
;     (see Pilyugin+06)
;   silent - do not print messages to STDOUT
;
; OUTPUTS:
;   abund - output data structure
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;    J. Moustakas, 2004 Feb 8, U of A, written
;
; Copyright (C) 2004, John Moustakas
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

function im_ff_abundance, line, snrcut_abundance=snrcut_abundance, nmonte=nmonte, $
  electrondensity=electrondensity, observed=observed, silent=silent
    
    nspec = n_elements(line)
    if (nspec eq 0L) then begin
       doc_library, 'im_ff_abundance'
       return, -1
    endif
    
    if (n_elements(snrcut_abundance) eq 0L) then snrcut_abundance = 3.0
    if not keyword_set(silent) then splog, 'S/N > '+string(snrcut_abundance,format='(G0.0)')+'.'

    if (n_elements(nmonte) eq 0L) then nmonte = 500
    if (n_elements(electrondensity) eq 0L) then electrondensity = 0

; initialize the output data structure    
    abund = {$
      ff_r23:                       -999.0, $ ; ([O II] + [O III] 4959,5007) / H-beta
      ff_r23_err:                   -999.0, $
      ff_o32:                       -999.0, $ ; [O III] 4959,5007 / [O II] 3727
      ff_o32_err:                   -999.0, $
      ff_r:                         -999.0, $ ; observed [O III] 4363 / H-beta
      ff_r_err:                     -999.0, $
      ff_r3:                        -999.0, $ ; [O III] 4959,5007 / H-beta
      ff_r3_err:                    -999.0, $
      ff_r2:                        -999.0, $ ; [O II] 3727 / H-beta
      ff_r2_err:                    -999.0, $
      ff_P:                         -999.0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      ff_P_err:                     -999.0, $
      ff_sii:                       -999.0, $ ; [S II] 6716 / [S II] 6731
      ff_sii_err:                   -999.0, $ 

      ff_r_pred:                    -999.0, $ ; predicted [O III] 4363 / H-beta
      ff_r_pred_err:                -999.0, $
      ff_t3:                        -999.0, $ ; T(O++) in units of 10,000 K
      ff_t3_err:                    -999.0, $
      ff_t2:                        -999.0, $ ; T(O+) in units of 10,000 K
      ff_t2_err:                    -999.0, $
      ff_log12oh_pp:                -999.0, $ ; 12+log(O++/H+)
      ff_log12oh_pp_err:            -999.0, $
      ff_log12oh_p:                 -999.0, $ ; 12+log(O+/H+)
      ff_log12oh_p_err:             -999.0, $
      ff_log12oh:                   -999.0, $ ; total 12+log(O/H)
      ff_log12oh_err:               -999.0, $
      ff_density:                    100.0}   ; electron density (default 100 cm-3)
              
    abund = replicate(abund,nspec)

; ---------------------------------------------------------------------------
; compute observed line ratios
; ---------------------------------------------------------------------------    

; r23    
    
    if (tag_exist(line[0],'OII_3727') and tag_exist(line[0],'OIII_4959') and $
      tag_exist(line[0],'OIII_5007') and tag_exist(line[0],'H_BETA')) then begin

       lineratio, line, ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          abund[index].ff_r23 = x
          abund[index].ff_r23_err = xerr
       endif

    endif

; o32

    if (tag_exist(line[0],'OII_3727') and tag_exist(line[0],'OIII_4959') and $
      tag_exist(line[0],'OIII_5007')) then begin

       lineratio, line, ['OIII_4959','OIII_5007'], 'OII_3727', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
    
       if (nindex ne 0L) then begin
          abund[index].ff_o32 = x
          abund[index].ff_o32_err = xerr
       endif

    endif

; observed R=[OIII] 4363/H-beta
    if (tag_exist(line[0],'OIII_4363') and tag_exist(line[0],'H_BETA')) then begin

       lineratio, line, 'OIII_4363', 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog
    
       if (nindex ne 0L) then begin
          abund[index].ff_r = x
          abund[index].ff_r_err = xerr
       endif

    endif

; r3
    
    if (tag_exist(line[0],'OIII_4959') and tag_exist(line[0],'OIII_5007') and $
      tag_exist(line[0],'H_BETA')) then begin

       lineratio, line, ['OIII_4959','OIII_5007'], 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          abund[index].ff_r3 = x
          abund[index].ff_r3_err = xerr
       endif

    endif

; r2
    
    if (tag_exist(line[0],'OII_3727') and tag_exist(line[0],'H_BETA')) then begin

       lineratio, line, 'OII_3727', 'H_BETA', '', '', x, xerr, $
         index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          abund[index].ff_r2 = x
          abund[index].ff_r2_err = xerr
       endif

    endif

; P-parameter
    
    if (tag_exist(line[0],'OIII_4959') and tag_exist(line[0],'OIII_5007') and $
      tag_exist(line[0],'OII_3727')) then begin

       lineratio, line, ['OIII_4959','OIII_5007'], ['OII_3727','OIII_4959','OIII_5007'], $
         '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          abund[index].ff_P = x
          abund[index].ff_P_err = xerr
       endif

    endif

; [S II] 6716 / [S II] 6731
    
    if (tag_exist(line[0],'SII_6716') and tag_exist(line[0],'SII_6731')) then begin

       lineratio, line, 'SII_6716', 'SII_6731', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          abund[index].ff_sii = x
          abund[index].ff_sii_err = xerr
       endif

    endif

; ---------------------------------------------------------------------------
; electron density [cm-3]
; ---------------------------------------------------------------------------

    index = where((abund.ff_sii gt -900.0),nindex)
    if (nindex ne 0L) and electrondensity then begin

       if (not keyword_set(silent)) then splog, 'Computing electron densities.'
       fivel = idl_fivel(1,11,lineratio=abund[index].ff_sii,$;err_lineratio=abund[index].ff_sii_err,$
         temperature=temperature,density=density,/silent)
       abund[index].ff_density = fivel.density
;      abund[index].ff_density_err = djs_mean([fivel.density_lower,fivel.density_upper])

    endif

; ---------------------------------------------------------------------------
; predicted R (Pilyugin et al. 2006, eq. 14); estimate the uncertainty
; in R empirically; NB! only appropriate for 12+log(O/H)>8.25
; ---------------------------------------------------------------------------

    index = where((abund.ff_p gt -900.0) and (abund.ff_r3 gt -900.0),nindex)
    if (nindex ne 0L) then begin

       p = abund[index].ff_p
       p_err = abund[index].ff_p_err
       p_monte = rebin(reform(p,1,nindex),nmonte,nindex)+randomn(seed,nmonte,nindex)*$
         rebin(reform(p_err,1,nindex),nmonte,nindex)

       r3 = abund[index].ff_r3
       r3_err = abund[index].ff_r3_err
       r3_monte = rebin(reform(r3,1,nindex),nmonte,nindex)+randomn(seed,nmonte,nindex)*$
         rebin(reform(r3_err,1,nindex),nmonte,nindex)

       c = [-4.151D,-3.118D,+2.958D,-0.680D]
       
       logr = c[0] + c[1]*alog10(p) + c[2]*alog10(r3) + c[3]*(alog10(p))^2.0
;      logr = -4.264 + 3.087*abund[index].ff_r23
       logr_err = logr*0.0

       for ii = 0L, nindex-1L do begin
          pos = where((p_monte[*,ii] gt 0.0) and (r3_monte[*,ii] gt 0.0) and (p_monte[*,ii] gt 0.0))
          r_monte = c[0] + c[1]*alog10(P_monte[pos,ii]) + c[2]*$
            alog10(r3_monte[pos,ii]) + c[3]*(alog10(p_monte[pos,ii]))^2.0
          logr_err[ii] = stddev(r_monte)
       endfor

       abund[index].ff_r_pred = 10.0^logr
       abund[index].ff_r_pred_err = alog(10.0)*logr_err*abund[index].ff_r_pred

;      plot, alog10(line[index].oiii_4959[0]+line[index].oiii_5007[0]), alog10(line[index].oiii_4363[0]), $
;        ps=4, xr=[-0.9,1.1], yr=[-5,0.0], xsty=3, ysty=3

;      niceprint, line[index].oiii_4363[0]*100, abund[index].ff_r*100.0

    endif

; ---------------------------------------------------------------------------
; iteratively compute the electron temperatures, t3 and t2, and the
; oxygen abundances in the two zones; labeled equations are from
; Pilyugin et al. (2006); by default use the predicted value of R
; given by Pilyugin, otherwise use the observed values
; ---------------------------------------------------------------------------

    if keyword_set(observed) then begin
       rratio = abund.ff_r
       rratio_err = abund.ff_r_err
    endif else begin
       rratio = abund.ff_r_pred
       rratio_err = abund.ff_r_pred_err
    endelse
    
    index = where((rratio gt -900.0) and (abund.ff_r3 gt -900.0) and (abund.ff_r2 gt -900.0),nindex)
    if (nindex ne 0L) then begin

; loop on each object

       for i = 0L, nindex-1L do begin
          
          r = rratio[index[i]]
          r_err = rratio_err[index[i]]
          r_monte = r + randomn(seed,nmonte)*r_err

          r3 = abund[index[i]].ff_r3
          r3_err = abund[index[i]].ff_r3_err
          r3_monte = r3 + randomn(seed,nmonte)*r3_err

          r2 = abund[index[i]].ff_r2
          r2_err = abund[index[i]].ff_r2_err
          r2_monte = r2 + randomn(seed,nmonte)*r2_err

          nelectron = abund[index[i]].ff_density ; [cm^-3]
       
          t3 = 1.0      ; initial guess (10,000 K)
          t3_old = 2.0

          while (abs(t3-t3_old)/t3 gt 0.001) do begin ; iterate

             t2 = 0.7*t3 + 0.3     ; eq. (9) [also, Campbell et al. (1986)]
      
             x2 = 1D-4*nelectron*t2^(-0.5)                   ; eq. (7)
             x3 = 1D-4*nelectron*t3^(-0.5)                   ; eq. (4)
             v = (1.0D + 0.0004D*x3) / (1.0D + 0.044D*x3)       ; eq. (3)
             CT = (8.44D - 1.09D*t3 + 0.5D*t3^2 - 0.08D*t3^3)*v ; eq. (2)
      
             t3_old = t3
             t3 = 1.432D / (alog10(r3/r) - alog10(CT)) ; eq (1)
;            print, t3, t3_old, abs(t3-t3_old)/t3 

; on the last iterature compute uncertainties; remove negative R3 and
; R values as being unphysical
             
             if (abs(t3-t3_old)/t3 le 0.001) then begin
                pos = where((r3_monte gt 0.0) and (r_monte gt 0.0))
                t3_err = stddev(1.432D / (alog10(r3_monte[pos]/r_monte[pos]) - alog10(CT)))
                t2_err = 0.7*t3_err
             endif
             
          endwhile

; compute the oxygen abundances and uncertainties
          
          t3_monte = t3 + randomn(seed,nmonte)*t3_err
          t2_monte = t2 + randomn(seed,nmonte)*t2_err
          
          log12oh_pp = alog10(r3) + 6.200D + 1.251D/t3 - 0.55D*alog10(t3) - 0.014D*t3                       ; eq (5)
          log12oh_p = alog10(r2) + 5.961D + 1.676D/t2 - 0.40D*alog10(t2) - 0.034D*t2 + alog10(1.0+1.35D*x2) ; eq (6)

          oh_pp = 10.0^(log12oh_pp - 12.0)
          oh_p = 10.0^(log12oh_p - 12.0)

          oh = oh_pp + oh_p ; eq (8)
          log12oh = alog10(oh) + 12.0

          pos = where((r3_monte gt 0.0) and (t3_monte gt 0.0))
          log12oh_pp_err = stddev(alog10(r3_monte[pos]) + 6.200D + 1.251D/t3_monte[pos] - $
            0.55D*alog10(t3_monte[pos]) - 0.014D*t3_monte[pos])
          pos = where((r2_monte gt 0.0) and (t2_monte gt 0.0))
          log12oh_p_err = stddev(alog10(r2_monte[pos]) + 5.961D + 1.676D/t2_monte[pos] - $
            0.40D*alog10(t2_monte[pos]) - 0.034D*t2_monte[pos] + alog10(1.0+1.35D*x2))

          oh_pp_err = alog(10.0)*log12oh_pp_err*oh_pp
          oh_p_err = alog(10.0)*log12oh_p_err*oh_p

          oh_err = sqrt(oh_pp_err^2.0 + oh_p_err^2.0)
          log12oh_err = oh_err/oh/alog(10.0)
          
          abund[index[i]].ff_t2 = t2
          abund[index[i]].ff_t3 = t3
          abund[index[i]].ff_log12oh_pp = log12oh_pp
          abund[index[i]].ff_log12oh_p = log12oh_p
          abund[index[i]].ff_log12oh = log12oh

          abund[index[i]].ff_t2_err = t2_err
          abund[index[i]].ff_t3_err = t3_err
          abund[index[i]].ff_log12oh_pp_err = log12oh_pp_err
          abund[index[i]].ff_log12oh_p_err = log12oh_p_err
          abund[index[i]].ff_log12oh_err = log12oh_err

;         help, abund[index[i]], /str & cc = get_kbrd(1)
       
       endfor

    endif

;   struct_print, struct_trimtags(abund,select='*log12oh*')

return, abund
end
