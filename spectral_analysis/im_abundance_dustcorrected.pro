;+
; NAME:
;       IM_ABUNDANCE()
;
; PURPOSE:
;       Compute empirical abundances and various relevant line-ratios
;       from strong emission lines. 
;
; CALLING SEQUENCE:
;       abund = im_abundance(line,snrcut_abundance=,/normalize,/HbHg,$
;          /combination,/electrondensity,/silent)
;
; INPUTS:
;       line   - input structure of de-reddened line fluxes
;
; OPTIONAL INPUTS:
;       snrcut_abundance - require that individual emission lines have
;                          a signal-to-noise ratio greater than
;                          SNRCUT_ABUNDANCE 
;
; KEYWORD PARAMETERS:
;       normalize       - divide all the emission-line fluxes and flux  
;                         errors by the flux at H-beta
;       HbHg            - by default, a well-measured EBV_HAHB is
;                         required to compute abundances for a given
;                         object; setting this keyword looks for a
;                         well-measured EBV_HBHG
;       combination     - *either* EBV_HAHB *or* EBV_HBHG must be
;                         well-measured
;       electrondensity - compute the electron density using the
;                         sulfur line-ratio
;       silent          - do not print messages to STDOUT
;
; OUTPUTS:
;       abund - output data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       LINERATIO, IDL_FIVEL(), DJS_MEAN(), IM_PREDICT_TOII()
;
; COMMENTS:
;       Be sure to first use ICLASSIFICATION() to remove AGN from the
;       sample.  This routine must be passed *de-reddened* line fluxes
;       from IUNRED_LINEDUST().
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;        J. Moustakas, 2004 Feb 8, U of A, written
;        jm04jul06uofa - added empirical Te-based abundances 
;        jm04jul22uofa - added SILENT and COMBINATION keywords 
;        jm04jul27uofa - added a warning flag if R23>1.1; added
;                        ELECTRONDENSITY keyword
;        jm04nov03uofa - documentation updated; HBETA keyword replaced
;                        with HbHg keyword for consistency with
;                        IUNRED_LINEDUST() 
;        jm04dec08uofa - general development of various empirical
;                        diagnostics; remove objects with log R23 > 1 
;        jm05jan01uofa - compute EW line-ratios and abundances
;        jm05jan06uofa - do not require reddening-corrected line
;                        fluxes; 
;-

function im_abundance, line, snrcut_abundance=snrcut_abundance, $
  normalize=normalize, HbHg=HbHg, combination=combination, $
  electrondensity=electrondensity, silent=silent
    
    light = 2.99792458D10 ; speed of light [cm/s]

    oratio = 2.984
    ocor = 1.0+1.0/oratio

    nspec = n_elements(line)
    if (nspec eq 0L) then begin
       print, 'Syntax - abund = im_abundance(line,snrcut_abundance=,/normalize,$'
       print, '   /HbHg,/combination,/electrondensity,/silent)'
       return, -1L
    endif
    
    if n_elements(snrcut_abundance) eq 0L then snrcut_abundance = 5.0
    if not keyword_set(silent) then splog, 'S/N > '+string(snrcut_abundance,format='(G0.0)')+'.'

    if n_elements(electrondensity) eq 0L then electrondensity = 0L

; initialize the output data structure    
    
    abund = {$
      Z_R23:                       -999D0, $ ; ([O II] + [O III] 4959,5007) / H-beta
      Z_R23_err:                   -999D0, $
      Z_O32:                       -999D0, $ ; [O III] 4959,5007 / [O II] 3727
      Z_O32_err:                   -999D0, $
      Z_O3N2:                      -999D0, $ ; ([O III] 5007 / H-beta) / ([N II] 6584 / H-alpha)
      Z_O3N2_err:                  -999D0, $
      Z_R3:                        -999D0, $ ; [O III] 4959,5007 / H-beta
      Z_R3_err:                    -999D0, $
      Z_P:                         -999D0, $ ; [O III] 4959,5007 / ([O II] + [O III] 4959,5007) 
      Z_P_err:                     -999D0, $
      Z_sii:                       -999D0, $ ; [S II] 6716 / [S II] 6731
      Z_sii_err:                   -999D0, $ 
      Z_nii_ha:                    -999D0, $ ; [N II] 6584 / H-alpha
      Z_nii_ha_err:                -999D0, $
      Z_nii_oii:                   -999D0, $ ; [N II] 6584 / [O II] 3727
      Z_nii_oii_err:               -999D0, $
      Z_nii_sii:                   -999D0, $ ; [N II] 6584 / [S II] 6716,6731
      Z_nii_sii_err:               -999D0, $
      Z_oii_oiii:                  -999D0, $ ; [O II] 3727 / [O III] 5007
      Z_oii_oiii_err:              -999D0, $
      Z_oiii_h_beta:               -999D0, $ ; [O III] 5007 / H-beta
      Z_oiii_h_beta_err:           -999D0, $
      Z_density:                     100D, $ ; electron density (default 100 cm-3)
;     Z_density_err:               -999D0, $
;     Z_T_oiii:                    -999D0, $ ; empirical T(O++) based on [N II]/Ha
;     Z_T_oiii_err:                -999D0, $
;     Z_T_oii:                     -999D0, $ ; empirical T(O+) from T(O++) and photoionization models
;     Z_T_oii_err:                 -999D0, $
;     Z_log12oh:                   -999D0, $ ; empirical 12+log (O/H) based on empirical T(O++)
;     Z_log12oh_err:               -999D0, $
;     Z_log_no:                    -999D0, $ ; empirical log (N/O) based on empirical T(O++)
;     Z_log_no_err:                -999D0, $
      Z_12OH_N2:                   -999D0, $ ; empirical abundance based my [N II]/Ha - Te relation
      Z_12OH_N2_err:               -999D0, $ 
      Z_12OH_N2_melbourne:         -999D0, $ ; Melbourne & Salzer (2002) [N II]/Ha calibration
      Z_12OH_N2_melbourne_err:     -999D0, $ 
      Z_12OH_N2_denicolo:          -999D0, $ ; Denicolo (2002) [N II]/Ha calibration
      Z_12OH_N2_denicolo_err:      -999D0, $ 
      Z_12OH_N2_pettini:           -999D0, $ ; Pettini & Pagel (2004) [N II]/Ha calibration
      Z_12OH_N2_pettini_err:       -999D0, $ 
      Z_12OH_O32:                  -999D0, $ ; empirical abundance based my O32 - Te relation
      Z_12OH_O32_err:              -999D0, $ 
      Z_12OH_P:                    -999D0, $ ; empirical abundance based my P - Te relation
      Z_12OH_P_err:                -999D0, $ 
      Z_12OH_O3N2:                 -999D0, $ ; empirical abundance based my O3N2 - Te relation - WRONG!!!
      Z_12OH_O3N2_err:             -999D0, $ 
      Z_12OH_O3N2_pettini:         -999D0, $ ; empirical abundance based on the Pettini & Pagel (2004) calibration
      Z_12OH_O3N2_pettini_err:     -999D0, $ 
      Z_12OH_CL01A:                -999D0, $ ; 12 + log O/H, Case A
      Z_12OH_CL01A_err:            -999D0, $
      Z_12OH_CL01F:                -999D0, $ ; 12 + log O/H, Case F
      Z_12OH_CL01F_err:            -999D0, $
      Z_12OH_ZKH94:                -999D0, $ ; R23 - upper branch only 
      Z_12OH_ZKH94_err:            -999D0, $
      Z_12OH_M91_upper:            -999D0, $ ; R23 - upper branch
      Z_12OH_M91_upper_err:        -999D0, $
      Z_12OH_M91_lower:            -999D0, $ ; R23 - lower branch
      Z_12OH_M91_lower_err:        -999D0, $
      Z_12OH_M91_contini:          -999D0, $ ; Contini et al. (2002) procedure for M91
      Z_12OH_M91_contini_err:      -999D0, $
      Z_12OH_M91_KD02:             -999D0, $ ; Kewley & Dopita (2002) procedure for M91
      Z_12OH_M91_KD02_err:         -999D0, $
      Z_12OH_M91_O32:              -999D0, $ ; use O32 to choose the branch
      Z_12OH_M91_O32_err:          -999D0, $
      Z_12OH_T04:                  -999D0, $ ; Tremonti et al. (2004)
      Z_12OH_T04_err:              -999D0, $
      Z_12OH_KK04_R23_upper:       -999D0, $ ; Kobulnicky & Kewley (2004) R23 upper branch
      Z_12OH_KK04_R23_upper_err:   -999D0, $
      Z_12OH_KK04_R23_lower:       -999D0, $ ; Kobulnicky & Kewley (2004) R23 lower branch
      Z_12OH_KK04_R23_lower_err:   -999D0, $
      Z_12OH_KK04_logU_upper:      -999D0, $ ; Kobulnicky & Kewley (2004) ionization parameter
      Z_12OH_KK04_logU_lower:      -999D0, $
      Z_12OH_KD02_nii_oii:         -999D0, $ ; Kewley & Dopita (2002) using [N II]/[O II]
      Z_12OH_KD02_nii_oii_err:     -999D0, $
;     Z_12OH_KD02_R23:             -999D0, $ ; Kewley & Dopita (2002) using R23
;     Z_12OH_KD02_R23_err:         -999D0, $
      Z_12OH_KD02_combined:        -999D0, $ ; Kewley & Dopita (2002) combined abundance diagnostic
      Z_12OH_KD02_combined_err:    -999D0, $
      Z_12OH_P01_upper:            -999D0, $ ; Pilyugin 2001, A&A, 369, 594
      Z_12OH_P01_upper_err:        -999D0, $
      Z_12OH_P01_lower:            -999D0, $
      Z_12OH_P01_lower_err:        -999D0, $
      Z_12OH_P01_average:          -999D0, $ ; average of the lower and upper calibrations, to be used in the turn-around region
      Z_12OH_P01_average_err:      -999D0, $
      Z_12OH_P01_melbourne:        -999D0, $ ; Melbourne & Salzer (2002) procedure for P01
      Z_12OH_P01_melbourne_err:    -999D0, $
      Z_12OH_P01_N2:               -999D0, $ ; use N2 to select the appropriate P01 calibration, using the average 
      Z_12OH_P01_N2_err:           -999D0, $ ;    of the two calibrations between 7.95 and 8.2
      Z_12OH_O3N2_N2:              -999D0, $ ; use O3N2 and N2 [Pettini] above and below 12+log(O/H)=8.5, respectively
      Z_12OH_O3N2_N2_err:          -999D0, $
      Z_EW_R23:                    -999D0, $ ; {EW([O II]) + EW([O III] 4959,5007)} / EW(H-beta)
      Z_EW_R23_err:                -999D0, $
      Z_EW_O32:                    -999D0, $ ; EW([O III] 4959,5007) / EW([O II] 3727)
      Z_EW_O32_err:                -999D0, $
      Z_EW_nii_ha:                 -999D0, $ ; EW([N II] 6584) / EW(H-alpha)
      Z_EW_nii_ha_err:             -999D0, $
      Z_EW_12OH_M91_upper:         -999D0, $ ; EW(R23) - upper branch
      Z_EW_12OH_M91_upper_err:     -999D0, $
      Z_EW_12OH_M91_lower:         -999D0, $ ; EW(R23) - lower branch
      Z_EW_12OH_M91_lower_err:     -999D0, $
      Z_EW_12OH_M91_O32:           -999D0, $ ; use EW(O32) to choose the branch
      Z_EW_12OH_M91_O32_err:       -999D0, $
      Z_EW_12OH_KK04_R23_upper:    -999D0, $ ; Kobulnicky & Kewley (2004) R23 upper branch
      Z_EW_12OH_KK04_R23_upper_err:-999D0, $
      Z_EW_12OH_KK04_R23_lower:    -999D0, $ ; Kobulnicky & Kewley (2004) R23 lower branch
      Z_EW_12OH_KK04_R23_lower_err:-999D0, $
      Z_EW_12OH_KK04_logU_upper:   -999D0, $ ; Kobulnicky & Kewley (2004) ionization parameter
      Z_EW_12OH_KK04_logU_lower:   -999D0}
        
    abund = replicate(abund,nspec)

; ###########################################################################
; compute EW line-ratios and abundances
; ###########################################################################

; EW(R23)

    if (tag_exist(line[0],'OII_3727_EW') and tag_exist(line[0],'OIII_4959_EW') and $
      tag_exist(line[0],'OIII_5007_EW') and tag_exist(line[0],'H_BETA_EW')) then begin

       lineratio, line, ['OII_3727_EW','OIII_4959_EW','OIII_5007_EW'], 'H_BETA_EW', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin

          abund[index].Z_EW_R23 = x
          abund[index].Z_EW_R23_err = xerr

;         toobig = where(x gt 1.0,ntoobig,comp=perfect,ncomp=nperfect)
;         if (ntoobig ne 0L) then begin
;            splog, 'Removing '+string(ntoobig,format='(I0.0)')+' objects '+$
;              'with log EW(R23) > 1.0.'
;         endif
;         if (nperfect ne 0L) then begin
;            abund[index[perfect]].Z_EW_R23 = x[perfect]
;            abund[index[perfect]].Z_EW_R23_err = xerr[perfect]
;         endif

       endif

    endif

; EW(O32)

    if (tag_exist(line[0],'OII_3727_EW') and tag_exist(line[0],'OIII_4959_EW') and $
      tag_exist(line[0],'OIII_5007_EW')) then begin

       lineratio, line, ['OIII_4959_EW','OIII_5007_EW'], 'OII_3727_EW', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
    
       if (nindex ne 0L) then begin
          abund[index].Z_EW_O32 = x
          abund[index].Z_EW_O32_err = xerr
       endif

    endif

    if (tag_exist(line[0],'NII_6584_EW') and tag_exist(line[0],'H_ALPHA_EW')) then begin

       lineratio, line, 'NII_6584_EW', 'H_ALPHA_EW', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          abund[index].Z_EW_nii_ha = x
          abund[index].Z_EW_nii_ha_err = xerr
       endif

    endif
    
; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999) - based on EW's
; ---------------------------------------------------------------------------    

    index = where((abund.Z_EW_R23 gt -900.0) and (abund.Z_EW_O32 gt -900.0),nindex)
    if (nindex ne 0L) then begin

; assume fixed uncertainties in the calibrations       
       
       x = abund[index].Z_EW_R23
       xerr = abund[index].Z_EW_R23_err

       y = abund[index].Z_EW_O32
       yerr = abund[index].Z_EW_O32_err

; lower branch    
    
       abund[index].Z_EW_12OH_M91_lower = 12.0 - 4.944 + 0.767*x + $
         0.602*x^2 - y*(0.29 + 0.332*x - 0.332*x^2)
       abund[index].Z_EW_12OH_M91_lower_err = 0.1

; force the lower-branch abundances to be less than 8.4

;      toobig = where(abund[index].Z_EW_12OH_M91_lower gt 8.4,ntoobig)
;      if (ntoobig ne 0L) then begin
;         abund[index[toobig]].Z_EW_12OH_M91_lower = -999.0
;         abund[index[toobig]].Z_EW_12OH_M91_lower_err = -999.0
;      endif
       
; upper branch    

       abund[index].Z_EW_12OH_M91_upper = 12.0 - 2.939 - 0.2*x - $
         0.237*x^2 - 0.305*x^3 - 0.0283*x^4 - $
         y*(0.0047 - 0.0221*x - 0.102*x^2 - 0.0817*x^3 - $
         0.00717*x^4)
       abund[index].Z_EW_12OH_M91_upper_err = 0.15
       
; force the upper-branch abundances to be greater than 8.4

;      toosmall = where(abund[index].Z_EW_12OH_M91_upper lt 8.4,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         abund[index[toosmall]].Z_EW_12OH_M91_upper = -999.0
;         abund[index[toosmall]].Z_EW_12OH_M91_upper_err = -999.0
;      endif
       
    endif

; ---------------------------------------------------------------------------    
; use EW(O32) to select the appropriate branch of the M91 calibration;
; in the turn-around region use my empirical O32 calibration (based on
; the line-fluxes); define the turn-around region as being between 8.2
; and 8.5 (based on the O32-12+log(O/H) calibration
; ---------------------------------------------------------------------------    

    index = where((abund.Z_EW_O32 gt -900.0) and (abund.Z_EW_12OH_M91_upper gt -900.0) and $
      (abund.Z_EW_12OH_M91_lower gt -900.0),nindex)
    if (nindex ne 0L) then begin

       lower = where(abund[index].Z_EW_O32 gt 0.6,nlower)
       if (nlower ne 0L) then begin
          abund[index[lower]].Z_EW_12OH_M91_O32 = abund[index[lower]].Z_EW_12OH_M91_lower
          abund[index[lower]].Z_EW_12OH_M91_O32_err = abund[index[lower]].Z_EW_12OH_M91_lower_err
       endif

       upper = where(abund[index].Z_EW_O32 lt 0.1,nupper)
       if (nupper ne 0L) then begin
          abund[index[upper]].Z_EW_12OH_M91_O32 = abund[index[upper]].Z_EW_12OH_M91_upper
          abund[index[upper]].Z_EW_12OH_M91_O32_err = abund[index[upper]].Z_EW_12OH_M91_upper_err
       endif

       turn = where((abund[index].Z_EW_O32 ge 0.1) and (abund[index].Z_EW_O32 le 0.6),nturn)
       if (nturn ne 0L) then begin
       
          c = [8.59D,-0.622D]
          cerr = [0.005D,0.010D]

          x = abund[index[turn]].Z_EW_O32
          xerr = abund[index[turn]].Z_EW_O32_err
          
          abund[index[turn]].Z_EW_12OH_M91_O32 = c[0] + c[1] * x
          abund[index[turn]].Z_EW_12OH_M91_O32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

       endif

    endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - R23
; ---------------------------------------------------------------------------    

    index = where((abund.Z_EW_O32 gt -900.0) and (abund.Z_EW_R23 gt -900.0),nindex)
    if (nindex ne 0L) then begin

       y = abund[index].Z_EW_O32
       x = abund[index].Z_EW_R23
       
; upper branch       

       logOH_upper = 9.0
       
       for j = 0L, 9L do begin ; iterate

          numer = 32.81D - 1.153*y^2 + logOH_upper*(-3.396 - 0.025*y + 0.1444*y^2)
          denom = 4.603D - 0.3119*y - 0.163*y^2 + logOH_upper*(-0.48 + 0.0271*y + 0.02037*y^2)
          logq = numer/denom

          logOH_upper = 9.72D - 0.777*x - 0.951*x^2 - 0.072*x^3 - 0.811*x^4 - $
            logq*(0.0737 - 0.0713*x - 0.141*x^2 + 0.0373*x^3 - 0.058*x^4)

       endfor

       abund[index].Z_EW_12OH_KK04_R23_upper = logOH_upper
       abund[index].Z_EW_12OH_KK04_logU_upper = logq - alog10(light)

; force the upper-branch abundances to be greater than 8.4

;      toosmall = where(abund[index].Z_EW_12OH_KK04_R23_upper lt 8.4,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         abund[index[toosmall]].Z_EW_12OH_KK04_R23_upper = -999.0
;         abund[index[toosmall]].Z_EW_12OH_KK04_R23_upper_err = -999.0
;         abund[index[toosmall]].Z_EW_12OH_KK04_logU_upper = -999.0
;      endif
       
; lower branch       
       
       logOH_lower = 7.9
       
       for j = 0L, 9L do begin ; iterate

          numer = 32.81D - 1.153*y^2 + logOH_lower*(-3.396 - 0.025*y + 0.1444*y^2)
          denom = 4.603D - 0.3119*y - 0.163*y^2 + logOH_lower*(-0.48 + 0.0271*y + 0.02037*y^2)
          logq = numer/denom

          logOH_lower = 9.40D + 4.65*x - 3.17*x^2 - logq*(0.272 + 0.547*x - 0.513*x^2)

       endfor

       abund[index].Z_EW_12OH_KK04_R23_lower = logOH_lower
       abund[index].Z_EW_12OH_KK04_logU_lower = logq - alog10(light)
       
; force the lower-branch abundances to be less than 8.4

;      toobig = where(abund[index].Z_EW_12OH_KK04_R23_lower gt 8.4,ntoobig)
;      if (ntoobig ne 0L) then begin
;         abund[index[toobig]].Z_EW_12OH_KK04_R23_lower = -999.0
;         abund[index[toobig]].Z_EW_12OH_KK04_R23_lower_err = -999.0
;         abund[index[toobig]].Z_EW_12OH_KK04_logU_lower = -999.0
;      endif
       
    endif

; ###########################################################################
; compute line-ratios and abundances
; ###########################################################################
    
    if keyword_set(HbHg) then begin

       if (tag_exist(line[0],'EBV_HBHG_ERR') eq 0L) then begin
          splog, 'Fluxes must be de-reddened.'
          return, abund
       endif
       
       good_ebv = where(line.ebv_hbhg_err gt 0.0,ngood)
       if (ngood eq 0L) then begin
          splog, 'No good de-reddened line fluxes.'
          return, abund
       endif

    endif 

    if keyword_set(combination) then begin

       if (tag_exist(line[0],'EBV_HAHB_ERR') eq 0L) and $
         (tag_exist(line[0],'EBV_HBHG_ERR') eq 0L) then begin
          splog, 'Fluxes must be de-reddened.'
          return, abund
       endif
       
       good_ebv = where((line.ebv_hahb_err gt 0.0) or (line.ebv_hbhg_err gt 0.0),ngood)
       if (ngood eq 0L) then begin
          splog, 'No good de-reddened line fluxes.'
          return, abund
       endif

    endif 

    if (not keyword_set(HbHg)) and (not keyword_set(combination)) then begin
       
       if (tag_exist(line[0],'EBV_HAHB_ERR') eq 0L) then begin
          splog, 'Fluxes must be de-reddened.'
          return, abund
       endif
       
       good_ebv = where(line.ebv_hahb_err gt 0.0,ngood)
       if (ngood eq 0L) then begin
          splog, 'No good de-reddened line fluxes.'
          return, abund
       endif

    endif
    
; if we need to normalize to H-beta then if H-beta is not
; well-measured then set that line to unmeasured as well; this should
; not be a big deal because we also cut on reddening, which requires
; H-beta; do not probagate the error in the normalization because that
; is done below
    
    if keyword_set(normalize) then begin

       tline = line
       
       linename = strtrim(tline[0].linename,2)
       nline = n_elements(linename)
       tags = tag_names(tline[0])

       hbtag = where(strmatch(tags,'H_BETA',/fold) eq 1B,nmatch)
       
       for j = 0L, nspec-1L do begin

          hbflux = reform((tline[j].(hbtag))[0,*])
          hbferr = reform((tline[j].(hbtag))[1,*])

          for k = 0L, nline-1L do begin 

             match = where(strmatch(tags,linename[k],/fold) eq 1B,nmatch)
             flux = (tline[j].(match))[0]
             ferr = (tline[j].(match))[1]
             if (ferr gt 0.0) and (hbferr gt 0.0) then begin
                flux = flux / hbflux
                ferr = ferr / hbflux
                tline[j].(match) = transpose([ [flux], [ferr] ])
             endif
             
          endfor

       endfor

       iline = tline[good_ebv]
    
    endif else iline = line[good_ebv]

    iabund = abund[good_ebv]

; ---------------------------------------------------------------------------
; compute relevant line ratios
; ---------------------------------------------------------------------------    

; R23    
    
    if (tag_exist(iline[0],'OII_3727') and tag_exist(iline[0],'OIII_4959') and $
      tag_exist(iline[0],'OIII_5007') and tag_exist(iline[0],'H_BETA')) then begin

       lineratio, iline, ['OII_3727','OIII_4959','OIII_5007'], 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          toobig = where(x gt 1.0,ntoobig,comp=perfect,ncomp=nperfect)
          if (ntoobig ne 0L) then begin
             splog, 'Removing '+string(ntoobig,format='(I0.0)')+' objects '+$
               'with log R23 > 1.0.'
          endif
          if (nperfect ne 0L) then begin
             iabund[index[perfect]].Z_R23 = x[perfect]
             iabund[index[perfect]].Z_R23_err = xerr[perfect]
          endif
       endif

    endif

; O32

    if (tag_exist(iline[0],'OII_3727') and tag_exist(iline[0],'OIII_4959') and $
      tag_exist(iline[0],'OIII_5007')) then begin

       lineratio, iline, ['OIII_4959','OIII_5007'], 'OII_3727', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance
    
       if (nindex ne 0L) then begin
          iabund[index].Z_O32 = x
          iabund[index].Z_O32_err = xerr
       endif

    endif

; O3N2 - Pettini & Pagel (2004)

    if (tag_exist(iline[0],'NII_6584') and tag_exist(iline[0],'OIII_5007') and $
      tag_exist(iline[0],'H_BETA') and tag_exist(iline[0],'H_ALPHA')) then begin

       lineratio, iline, 'OIII_5007', 'H_BETA', 'NII_6584', 'H_ALPHA', $
         x1, x1err, x2, x2err, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_O3N2 = x1 - x2
          iabund[index].Z_O3N2_err = sqrt(x1err^2 + x2err^2)
       endif

    endif

; R3 (Pilyugin 2001)
    
    if (tag_exist(iline[0],'OIII_4959') and tag_exist(iline[0],'OIII_5007') and $
      tag_exist(iline[0],'H_BETA')) then begin

       lineratio, iline, ['OIII_4959','OIII_5007'], 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_R3 = x
          iabund[index].Z_R3_err = xerr
       endif

    endif

; P-parameter (Pilyugin 2001)
    
    if (tag_exist(iline[0],'OIII_4959') and tag_exist(iline[0],'OIII_5007') and $
      tag_exist(iline[0],'OII_3727')) then begin

       lineratio, iline, ['OIII_4959','OIII_5007'], ['OII_3727','OIII_4959','OIII_5007'], $
         '', '', x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          iabund[index].Z_P = x
          iabund[index].Z_P_err = xerr
       endif

    endif

; [S II] 6716 / [S II] 6731
    
    if (tag_exist(iline[0],'SII_6716') and tag_exist(iline[0],'SII_6731')) then begin

       lineratio, iline, 'SII_6716', 'SII_6731', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance, /nolog

       if (nindex ne 0L) then begin
          iabund[index].Z_sii = x
          iabund[index].Z_sii_err = xerr
       endif

    endif

; [N II]/Ha    
    
    if (tag_exist(iline[0],'NII_6584') and tag_exist(iline[0],'H_ALPHA')) then begin

       lineratio, iline, 'NII_6584', 'H_ALPHA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_nii_ha = x
          iabund[index].Z_nii_ha_err = xerr
       endif

    endif
    
; [N II]/[O II]
    
    if (tag_exist(iline[0],'NII_6584') and tag_exist(iline[0],'OII_3727')) then begin

       lineratio, iline, 'NII_6584', 'OII_3727', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_nii_oii = x
          iabund[index].Z_nii_oii_err = xerr
       endif

    endif

; [N II]/[S II]

    if (tag_exist(iline[0],'NII_6584') and tag_exist(iline[0],'SII_6716') and $
      tag_exist(iline[0],'SII_6731')) then begin

       lineratio, iline, 'NII_6584', ['SII_6716','SII_6731'], '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_nii_sii = x
          iabund[index].Z_nii_sii_err = xerr
       endif

    endif

; [O II]/[O III]
    
    if (tag_exist(iline[0],'OII_3727') and tag_exist(iline[0],'OIII_5007')) then begin

       lineratio, iline, 'OII_3727', 'OIII_5007', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_oii_oiii = x
          iabund[index].Z_oii_oiii_err = xerr
       endif

    endif

; [O III]/H-beta
    
    if (tag_exist(iline[0],'H_BETA') and tag_exist(iline[0],'OIII_5007')) then begin

       lineratio, iline, 'OIII_5007', 'H_BETA', '', '', $
         x, xerr, index=index, nindex=nindex, snrcut=snrcut_abundance

       if (nindex ne 0L) then begin
          iabund[index].Z_oiii_h_beta = x
          iabund[index].Z_oiii_h_beta_err = xerr
       endif

    endif

; ---------------------------------------------------------------------------    
; electron density [cm-3]
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_sii gt -900.0),nindex)
    if (nindex ne 0L) and electrondensity then begin

       if (not keyword_set(silent)) then splog, 'Computing electron densities.'
       fivel = idl_fivel(1,11,lineratio=iabund[index].Z_sii,$;err_lineratio=iabund[index].Z_sii_err,$
         temperature=temperature,density=density,/silent)
       iabund[index].Z_density = fivel.density
;      iabund[index].Z_density_err = djs_mean([fivel.density_lower,fivel.density_upper])

    endif
    
;;; ---------------------------------------------------------------------------    
;;; compute an empirical T(O++) based on [N II]/Ha then predict T(O+) 
;;; from the photionization models
;;; ---------------------------------------------------------------------------    
;;
;;    index = where(iabund.Z_nii_ha gt -900.0,nindex)
;;    if (nindex ne 0L) then begin
;;
;;       c = [2512D,-7578D]
;;       cerr = [359D,239D]
;;
;;       x = iabund[index].Z_nii_ha
;;       xerr = iabund[index].Z_nii_ha_err
;;
;;       inrange = where((x gt -2.5) and (x lt -0.3),ninrange)
;;       if (ninrange ne 0L) then begin
;;          
;;          iabund[index[inrange]].Z_T_oiii = c[0] + c[1] * x[inrange]
;;          iabund[index[inrange]].Z_T_oiii_err = sqrt( cerr[0]^2 + (c[1] * xerr[inrange])^2.0 + $
;;            (cerr[1] * x[inrange])^2.0 )
;;
;;          t_oii = im_predict_toii(iabund[index[inrange]].Z_T_oiii,toiii_lower=iabund[index[inrange]].Z_T_oiii_err,$
;;            toiii_upper=iabund[index[inrange]].Z_T_oiii_err)
;;
;;          iabund[index[inrange]].Z_T_oii     = t_oii.ZT_T_oii
;;          iabund[index[inrange]].Z_T_oii_err = t_oii.ZT_T_oii_err
;;          
;;          oh12Te_input = struct_addtags(iline[index[inrange]],im_struct_trimtags(iabund[index[inrange]],$
;;            select=['Z_DENSITY','Z_T_OIII','Z_T_OIII_ERR','Z_T_OII','Z_T_OII_ERR'],$
;;            newtags=['ZT_DENSITY','ZT_T_OIII','ZT_T_OIII_ERR','ZT_T_OII','ZT_T_OII_ERR']))
;;
;;          splog, 'Computing empirical electron temperature abundances.'
;;          oh12Te = im_compute_oh12te(oh12Te_input)
;;
;;          iabund[index[inrange]].Z_log12oh     = oh12Te.ZT_log12oh
;;          iabund[index[inrange]].Z_log12oh_err = oh12Te.ZT_log12oh_err
;;          iabund[index[inrange]].Z_log_no      = oh12Te.ZT_log_no
;;          iabund[index[inrange]].Z_log_no_err  = oh12Te.ZT_log_no_err
;;
;;       endif
;;          
;;    endif
;;
; ---------------------------------------------------------------------------    
; 12+log (O/H) based on my [N II]/Ha-Te calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_nii_ha gt -900),nindex)
;   index = where((iabund.Z_nii_ha ge -2.5) and (iabund.Z_nii_ha le -0.3),nindex)
    if (nindex ne 0L) then begin

       x = iabund[index].z_nii_ha
       xerr = iabund[index].z_nii_ha_err

       c = [9.153,0.782]
       cerr = [0.0485,0.0307]

;      iabund[index].Z_12oh_n2 = c[0] + c[1] * x + c[2] * x^2
;      iabund[index].Z_12oh_n2_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + $
;        (cerr[1]*x)^2.0 + (2*c[2]*x*xerr)^2.0 + (cerr[2]*x^2)^2.0 )

       iabund[index].Z_12oh_n2 = c[0] + c[1] * x
       iabund[index].Z_12oh_n2_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

    endif
    
; ---------------------------------------------------------------------------    
; Melbourne & Salzer (2002) [N II]/Ha calibration
; ---------------------------------------------------------------------------    

    index = where(iabund.Z_nii_ha gt -900.0,nindex)
    if (nindex ne 0L) then begin

       c = [9.26D,1.23D,0.204D]
;      cerr = []

       x = iabund[index].z_nii_ha
;      xerr = iabund[index].z_nii_ha_err

       iabund[index].Z_12oh_n2_melbourne = c[0] + c[1] * x + c[2] * x^2
       iabund[index].Z_12oh_n2_melbourne_err = 0.156

    endif
    
; ---------------------------------------------------------------------------    
; Denicolo et al. (2002) [N II]/Ha calibration
; ---------------------------------------------------------------------------    

    index = where(iabund.Z_nii_ha gt -900.0,nindex)
    if (nindex ne 0L) then begin

       c = [9.12D,0.73D]
       cerr = [0.05D,0.10D]
       
       x = iabund[index].z_nii_ha
       xerr = iabund[index].z_nii_ha_err
       
       iabund[index].Z_12oh_n2_denicolo = c[0] + c[1] * x
       iabund[index].Z_12oh_n2_denicolo_err = sqrt( cerr[0]^2 + (c[1] * xerr)^2.0 + (cerr[1] * x)^2.0 )

    endif
    
; ---------------------------------------------------------------------------    
; Pettini & Pagel (2004) [N II]/Ha calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_nii_ha ge -2.5) and (iabund.Z_nii_ha lt -0.3),nindx)
    if (nindex ne 0L) then begin

       c = [8.90D,0.57D]
       cerr = [0.18D,0.0D]
       
       x = iabund[index].z_nii_ha
       xerr = iabund[index].z_nii_ha_err
       
       iabund[index].Z_12oh_n2_pettini = c[0] + c[1] * x
       iabund[index].Z_12oh_n2_pettini_err = sqrt( cerr[0]^2 + (c[1] * xerr)^2.0 + (cerr[1] * x)^2.0 )

    endif
    
; ---------------------------------------------------------------------------    
; 12+log (O/H) based on my O32-Te calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_O32 ge -0.9) and (iabund.Z_O32 le 1.2),nindex)
    if (nindex ne 0L) then begin

       c = [8.2835434,-0.47751407]
       cerr = [0.024664146,0.042504591]

       x = iabund[index].Z_O32
       xerr = iabund[index].Z_O32_err

       iabund[index].Z_12oh_o32 = c[0] + c[1] * x
       iabund[index].Z_12oh_o32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

    endif
    
; ---------------------------------------------------------------------------    
; 12+log (O/H) based on my P-Te calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_P ge 0.1) and (iabund.Z_P le 1.0),nindex)
;   index = where((iabund.Z_P ge -900.0),nindx)
    if (nindex ne 0L) then begin

       c = [8.8139881,-1.0486419]
       cerr = [0.067280689,0.095832449]

       x = iabund[index].Z_P
       xerr = iabund[index].Z_P_err

       iabund[index].Z_12oh_P = c[0] + c[1] * x
       iabund[index].Z_12oh_P_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

    endif
    
;;; ---------------------------------------------------------------------------    
;;; 12+log (O/H) based on my O3N2-Te calibration - WRONG!!!
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((iabund.Z_O3N2 ge -0.7) and (iabund.Z_O3N2 le 1.9),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       c = [8.5638200,-0.21535851]
;;       cerr = [0.046192925,0.032765583]
;;
;;       x = iabund[index].Z_O3N2
;;       xerr = iabund[index].Z_O3N2_err
;;
;;       iabund[index].Z_12oh_o3n2 = c[0] + c[1] * x
;;       iabund[index].Z_12oh_o3n2_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )
;;
;;    endif
    
; ---------------------------------------------------------------------------    
; 12+log (O/H) based on the Pettini & Pagel (2004) O3N2-Te calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_O3N2 ge -1.0) and (iabund.Z_O3N2 le 1.9),nindex)
    if (nindex ne 0L) then begin

       c = [8.73D,-0.32D]
       cerr = [0.14D,0.0D]

       x = iabund[index].Z_O3N2
       xerr = iabund[index].Z_O3N2_err

       iabund[index].Z_12oh_o3n2_pettini = c[0] + c[1] * x
       iabund[index].Z_12oh_o3n2_pettini_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

    endif
    
; ---------------------------------------------------------------------------    
; Charlot & Longhetti (2001, Case A)    
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_nii_sii gt -900.0) and (iabund.Z_oii_oiii gt -900.0),nindex)
    if (nindex ne 0L) then begin

       x2 = 10^(iabund[index].Z_oii_oiii)/1.5
       x6 = 10^(iabund[index].Z_nii_sii)/0.85
       
       iabund[index].Z_12OH_CL01A = 12.0 + alog10(5.09D-4 * x2^0.17 * x6^1.17)
       iabund[index].Z_12OH_CL01A_err = 0.24

    endif
    
; ---------------------------------------------------------------------------    
; Charlot & Longhetti (2001, Case F)    
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_oiii_h_beta gt -900.0) and (iabund.Z_oii_oiii gt -900.0),nindex)

    if (nindex ne 0L) then begin

       oii_oiii = 10^iabund[index].Z_oii_oiii
       x2 = oii_oiii/1.5
       x3 = 10^iabund[index].Z_oiii_h_beta/2.0
       
       lo = where(oii_oiii lt 0.8,nlo,comp=hi,ncomp=nhi)

       if (nlo ne 0L) then begin
          iabund[index[lo]].Z_12OH_CL01F = 12 + alog10(3.78D-4 * x2[lo]^(0.17) * x3[lo]^(-0.44))
          iabund[index[lo]].Z_12OH_CL01F_err = 0.51
       endif
       
       if (nhi ne 0L) then begin
          iabund[index[hi]].Z_12OH_CL01F = 12 + alog10(3.96D-4 * x3[hi]^(-0.46))
          iabund[index[hi]].Z_12OH_CL01F_err = 0.31
       endif

    endif
    
; ---------------------------------------------------------------------------    
; Zaritsky et al. (1994)
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_R23 gt -900.0),nindex)
    if (nindex ne 0L) then begin

; do not use the formal errors in R23 because they are too small;
; assign a fixed uncertainty of 0.2 dex       
       
       x = iabund[index].Z_R23
       xerr = iabund[index].Z_R23_err

       c = [9.265D,-0.33,-0.202,-0.207,-0.333]
       
       iabund[index].Z_12OH_ZKH94 = c[0] + c[1]*x + c[2]*x^2 +c[3]*x^3 + c[4]*x^4
;      iabund[index].Z_12OH_ZKH94_err = sqrt( (c[1]*xerr)^2 + (2*c[2]*x*xerr)^2 + $
;        (3*c[3]*x^2*xerr)^2 + (4*c[4]*x^3*xerr)^2 )
       iabund[index].Z_12OH_ZKH94_err = 0.2
       
; force the ZKH abundances to be greater than 8.4

;      toosmall = where(iabund[index].Z_12OH_ZKH94 lt 8.4,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         iabund[index[toosmall]].Z_12OH_ZKH94 = -999.0
;         iabund[index[toosmall]].Z_12OH_ZKH94_err = -999.0
;      endif
       
    endif

; ---------------------------------------------------------------------------    
; McGaugh (1991) as fitted in Kobulnicky et al. (1999)
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_R23 gt -900.0) and (iabund.Z_O32 gt -900.0),nindex)
    if (nindex ne 0L) then begin

; assume fixed uncertainties in the calibrations       
       
       x = iabund[index].Z_R23
       xerr = iabund[index].Z_R23_err

       y = iabund[index].Z_O32
       yerr = iabund[index].Z_O32_err

; lower branch    
    
       iabund[index].Z_12OH_M91_lower = 12.0 - 4.944 + 0.767*x + $
         0.602*x^2 - y*(0.29 + 0.332*x - 0.332*x^2)
       iabund[index].Z_12OH_M91_lower_err = 0.1
       
; force the lower-branch abundances to be less than 8.4; see comments
; on the upper branch below 

;      toobig = where(iabund[index].Z_12OH_M91_lower gt 8.4,ntoobig)
;      if (ntoobig ne 0L) then begin
;         iabund[index[toobig]].Z_12OH_M91_lower = -999.0
;         iabund[index[toobig]].Z_12OH_M91_lower_err = -999.0
;      endif
       
; upper branch    

       iabund[index].Z_12OH_M91_upper = 12.0 - 2.939 - 0.2*x - $
         0.237*x^2 - 0.305*x^3 - 0.0283*x^4 - $
         y*(0.0047 - 0.0221*x - 0.102*x^2 - 0.0817*x^3 - $
         0.00717*x^4)
       iabund[index].Z_12OH_M91_upper_err = 0.15
       
; force the upper-branch abundances to be greater than 8.4

; jm05jan04uofa: i have decided to remove this condition because of
; the dependence of the turn-around region on the ionization
; parameter; even when Z_R23=1.0, M91_upper never dips below 8.11

;      toosmall = where(iabund[index].Z_12OH_M91_upper lt 8.4,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         iabund[index[toosmall]].Z_12OH_M91_upper = -999.0
;         iabund[index[toosmall]].Z_12OH_M91_upper_err = -999.0
;      endif
       
    endif

; ---------------------------------------------------------------------------    
; use the procedure in Contini et al. (2002), Section 4.1, to
; determine if an object is in the upper or lower branch of the M91
; relation
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_nii_ha gt -900.0) and (iabund.Z_nii_oii gt -900.0) and $
      (iabund.Z_12OH_M91_upper gt -900.0),nindex)
    if (nindex ne 0L) then begin

; [N II]/Ha and [N II]/[O II] both indicate the lower branch
       
       lower = where((iabund[index].Z_nii_ha lt -1.0) and (iabund[index].Z_nii_oii lt -1.05),nlower)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_M91_contini = iabund[index[lower]].Z_12OH_M91_lower
          iabund[index[lower]].Z_12OH_M91_contini_err = iabund[index[lower]].Z_12OH_M91_lower_err
       endif

; [N II]/Ha and [N II]/[O II] both indicate the upper branch

       upper = where((iabund[index].Z_nii_ha gt -1.0) and (iabund[index].Z_nii_oii gt -0.8),nupper)
       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_M91_contini = iabund[index[upper]].Z_12OH_M91_upper
          iabund[index[upper]].Z_12OH_M91_contini_err = iabund[index[upper]].Z_12OH_M91_upper_err
       endif

; [N II]/[O II] is in an ambigious region --> use [N II]/Ha

       ambig = where((iabund[index].Z_nii_oii gt -1.05) and (iabund[index].Z_nii_oii lt -0.8),nambig)

       if (nambig ne 0L) then begin

          lower = where((iabund[index[ambig]].Z_nii_ha lt -1.0),nlower,comp=upper,ncomp=nupper)

          if (nlower ne 0L) then begin
             iabund[index[ambig[lower]]].Z_12OH_M91_contini = iabund[index[ambig[lower]]].Z_12OH_M91_lower
             iabund[index[ambig[lower]]].Z_12OH_M91_contini_err = iabund[index[ambig[lower]]].Z_12OH_M91_lower_err
          endif

          if (nupper ne 0L) then begin
             iabund[index[ambig[upper]]].Z_12OH_M91_contini = iabund[index[ambig[upper]]].Z_12OH_M91_upper
             iabund[index[ambig[upper]]].Z_12OH_M91_contini_err = iabund[index[ambig[upper]]].Z_12OH_M91_upper_err
          endif

       endif

; [N II]/Ha indicates upper branch but [N II]/[O II] indicates lower
; branch --> take an average of the two metallicities
       
       veryambig = where((iabund[index].Z_nii_ha gt -1.0) and (iabund[index].Z_nii_oii lt -1.05),nambig)

       if (nambig ne 0L) then begin

          Z1 = iabund[index[veryambig]].Z_12OH_M91_lower
          Z1err = iabund[index[veryambig]].Z_12OH_M91_lower_err

          Z2 = iabund[index[veryambig]].Z_12OH_M91_upper
          Z2err = iabund[index[veryambig]].Z_12OH_M91_upper_err
          
          iabund[index[veryambig]].Z_12OH_M91_contini = (Z1 + Z2) / 2.0
          iabund[index[veryambig]].Z_12OH_M91_contini_err = sqrt(Z1err^2 + Z2err^2)

       endif
          
; [N II]/Ha indicates lower branch but [N II]/[O II] indicates upper
; branch --> take an average of the two metallicities (should rarely
; occur) 
       
       veryambig = where((iabund[index].Z_nii_ha lt -1.0) and (iabund[index].Z_nii_oii gt -0.8),nambig)

       if (nambig ne 0L) then begin

          Z1 = iabund[index[veryambig]].Z_12OH_M91_lower
          Z1err = iabund[index[veryambig]].Z_12OH_M91_lower_err

          Z2 = iabund[index[veryambig]].Z_12OH_M91_upper
          Z2err = iabund[index[veryambig]].Z_12OH_M91_upper_err
          
          iabund[index[veryambig]].Z_12OH_M91_contini = (Z1 + Z2) / 2.0
          iabund[index[veryambig]].Z_12OH_M91_contini_err = sqrt(Z1err^2 + Z2err^2)

       endif

    endif

; ---------------------------------------------------------------------------    
; Kobulnicky & Kewley (2004) - R23
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_O32 gt -900.0) and (iabund.Z_R23 gt -900.0),nindex)
    if (nindex ne 0L) then begin

       y = iabund[index].Z_O32
       x = iabund[index].Z_R23
       
; upper branch       

       logOH_upper = 9.0
       
       for j = 0L, 9L do begin ; iterate

          numer = 32.81D - 1.153*y^2 + logOH_upper*(-3.396 - 0.025*y + 0.1444*y^2)
          denom = 4.603D - 0.3119*y - 0.163*y^2 + logOH_upper*(-0.48 + 0.0271*y + 0.02037*y^2)
          logq = numer/denom

          logOH_upper = 9.72D - 0.777*x - 0.951*x^2 - 0.072*x^3 - 0.811*x^4 - $
            logq*(0.0737 - 0.0713*x - 0.141*x^2 + 0.0373*x^3 - 0.058*x^4)

       endfor

       iabund[index].Z_12OH_KK04_R23_upper = logOH_upper
       iabund[index].Z_12OH_KK04_logU_upper = logq - alog10(light)

; force the upper-branch abundances to be greater than 8.4

;      toosmall = where(iabund[index].Z_12OH_KK04_R23_upper lt 8.4,ntoosmall)
;      if (ntoosmall ne 0L) then begin
;         iabund[index[toosmall]].Z_12OH_KK04_R23_upper = -999.0
;         iabund[index[toosmall]].Z_12OH_KK04_R23_upper_err = -999.0
;         iabund[index[toosmall]].Z_12OH_KK04_logU_upper = -999.0
;      endif
       
; lower branch       
       
       logOH_lower = 7.9
       
       for j = 0L, 9L do begin ; iterate

          numer = 32.81D - 1.153*y^2 + logOH_lower*(-3.396 - 0.025*y + 0.1444*y^2)
          denom = 4.603D - 0.3119*y - 0.163*y^2 + logOH_lower*(-0.48 + 0.0271*y + 0.02037*y^2)
          logq = numer/denom

          logOH_lower = 9.40D + 4.65*x - 3.17*x^2 - logq*(0.272 + 0.547*x - 0.513*x^2)

       endfor

       iabund[index].Z_12OH_KK04_R23_lower = logOH_lower
       iabund[index].Z_12OH_KK04_logU_lower = logq - alog10(light)
       
; force the lower-branch abundances to be less than 8.4

;      toobig = where(iabund[index].Z_12OH_KK04_R23_lower gt 8.4,ntoobig)
;      if (ntoobig ne 0L) then begin
;         iabund[index[toobig]].Z_12OH_KK04_R23_lower = -999.0
;         iabund[index[toobig]].Z_12OH_KK04_R23_lower_err = -999.0
;         iabund[index[toobig]].Z_12OH_KK04_logU_lower = -999.0
;      endif
       
    endif

; ---------------------------------------------------------------------------    
; use ZKH94 or CL01A to select the appropriate branch according to
; Kewley & Dopita (2002)
; ---------------------------------------------------------------------------    

;   index = where((iabund.Z_12OH_CL01A gt -900.0) and (iabund.Z_12OH_M91_upper gt -900.0),nindex)
    index = where((iabund.Z_12OH_ZKH94 gt -900.0) and (iabund.Z_12OH_M91_upper gt -900.0),nindex)
    if (nindex ne 0L) then begin

;      lower = where(iabund[index].Z_12OH_CL01A le 8.4,nlower,comp=upper,ncomp=nupper)
       lower = where(iabund[index].Z_12OH_ZKH94 le 8.4,nlower,comp=upper,ncomp=nupper)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_M91_KD02 = iabund[index[lower]].Z_12OH_M91_lower
          iabund[index[lower]].Z_12OH_M91_KD02_err = iabund[index[lower]].Z_12OH_M91_lower_err
       endif

       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_M91_KD02 = iabund[index[upper]].Z_12OH_M91_upper
          iabund[index[upper]].Z_12OH_M91_KD02_err = iabund[index[upper]].Z_12OH_M91_upper_err
       endif
       
    endif
    
; ---------------------------------------------------------------------------    
; Kewley & Dopita (2002) [N II]/[O II], use the coefficients from
; Table 3 at the mid-point of the ionization parameter range, log U =
; -3.17 (q=2E7), which should give the same abundance as equation 5
; above 8.6. 
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_nii_oii gt -900.0),nindex)
    if (nindex ne 0L) then begin

       x = iabund[index].Z_nii_oii
       coeff = [1106.87,-532.154,96.3733,-7.81061,0.239282] # (fltarr(nindex)+1) ; q = 2E7
       coeff[0,*] = coeff[0,*] - x

; find the roots of this polynomial; the correct one must be real and
; between 7.5 and 9.4

       for i = 0L, nindex-1L do begin
          roots = fz_roots(coeff[*,i],/double)
          for j = 0L, 3L do begin
             if (imaginary(roots[j]) eq 0.0) then begin
                if ((abs(roots[j]) ge 7.5) and (abs(roots[j]) le 9.4)) then $
                  iabund[index[i]].Z_12OH_KD02_nii_oii = abs(roots[j])
             endif
          endfor
       endfor
                
    endif

; ---------------------------------------------------------------------------    
; Kewley & Dopita (2002) combined diagnostic; use KD02 [N II]/[O II]
; if [N II]/[O II] gives 8.6 < 12 + log (O/H); use the average of M91
; and ZKH94 if 8.5 < 12 + log (O/H) < 8.6; finally, use the average of
; CL01A and KK04-R23 if log(O/H)+12 < 8.5.
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_12oh_kd02_nii_oii gt -900.0) and (iabund.Z_12oh_m91_upper gt -900.0) and $
      (iabund.Z_12oh_zkh94 gt -900.0) and (iabund.Z_12oh_kk04_R23_lower gt -900.0) and $
      (iabund.Z_12oh_cl01a gt -900.0),nindex)

    if (nindex ne 0L) then begin

       M91Z94_average = (iabund[index].Z_12oh_m91_upper + iabund[index].Z_12oh_zkh94) / 2.0
       KD02C01_average = (iabund[index].Z_12oh_kk04_R23_lower + iabund[index].Z_12oh_cl01a) / 2.0

       lo = where(iabund[index].Z_12oh_kd02_nii_oii le 8.6,nlo,comp=up,ncomp=nup)

       if (nlo ne 0L) then begin

          upper = where(M91Z94_average[lo] ge 8.5,nupper,comp=lower,ncomp=nlower)
          if (nupper ne 0L) then iabund[index[lo[upper]]].Z_12oh_kd02_combined = M91Z94_average[lo[upper]]
          if (nlower ne 0L) then iabund[index[lo[lower]]].Z_12oh_kd02_combined = KD02C01_average[lo[lower]]
          
       endif

       if (nup ne 0L) then begin

          iabund[index[up]].Z_12oh_kd02_combined = iabund[index[up]].Z_12oh_kd02_nii_oii
          
       endif
                
    endif

; ---------------------------------------------------------------------------    
; use O32 to select the appropriate branch of the M91 calibration; in
; the turn-around region use an empirical O32 calibration; define the
; turn-around region as being between 8.2 and 8.5 (based on the
; O32-12+log(O/H) calibration
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_O32 gt -900.0),nindex)
    if (nindex ne 0L) then begin

       lower = where((iabund[index].Z_O32 gt 0.6) and (iabund.Z_12OH_M91_lower gt -900.0),nlower)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_M91_O32 = iabund[index[lower]].Z_12OH_M91_lower
          iabund[index[lower]].Z_12OH_M91_O32_err = iabund[index[lower]].Z_12OH_M91_lower_err
       endif

       upper = where((iabund[index].Z_O32 lt 0.1) and (iabund.Z_12OH_M91_upper gt -900.0),nupper)
       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_M91_O32 = iabund[index[upper]].Z_12OH_M91_upper
          iabund[index[upper]].Z_12OH_M91_O32_err = iabund[index[upper]].Z_12OH_M91_upper_err
       endif

       turn = where((iabund[index].Z_O32 ge 0.1) and (iabund[index].Z_O32 le 0.6),nturn)
       if (nturn ne 0L) then begin
       
          c = [8.59D,-0.622D]
          cerr = [0.005D,0.010D]

          x = iabund[index[turn]].Z_O32
          xerr = iabund[index[turn]].Z_O32_err
          
          iabund[index[turn]].Z_12OH_M91_O32 = c[0] + c[1] * x
          iabund[index[turn]].Z_12OH_M91_O32_err = sqrt( cerr[0]^2 + (c[1]*xerr)^2.0 + (cerr[1]*x)^2.0 )

       endif

    endif

; ---------------------------------------------------------------------------    
; Tremonti et al. (2004) - R23, upper branch only
; ---------------------------------------------------------------------------    

    index = where(iabund.Z_R23 gt -900.0,nindex)
    if (nindex ne 0L) then begin

       x = iabund[index].Z_R23

       OH12 = 9.185D - 0.313*x - 0.264*x^2 - 0.321*x^3
       upper = where((OH12 ge 8.5),nupper) ; applicable range of the polynomial

       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_T04 = OH12[upper]
          iabund[index[upper]].Z_12OH_T04_err = 0.12
       endif

    endif
    
;;; ---------------------------------------------------------------------------    
;;; use [O III]/Hb to select the appropriate branch of the KK04 calibration
;;; ---------------------------------------------------------------------------    
;;
;;    index = where((iabund.Z_oiii_h_beta gt -900.0) and (iabund.Z_12OH_KK04_R23_upper gt -900.0) and $
;;      (iabund.Z_12OH_KK04_R23_lower gt -900.0),nindex)
;;    if (nindex ne 0L) then begin
;;
;;       lower = where(iabund[index].Z_oiii_h_beta gt 0.25,nlower,comp=upper,ncomp=nupper)
;;       if (nlower ne 0L) then iabund[index[lower]].Z_12OH_KK04_O32 = iabund[index[lower]].Z_12OH_KK04_R23_lower
;;       if (nupper ne 0L) then iabund[index[upper]].Z_12OH_KK04_O32 = iabund[index[upper]].Z_12OH_KK04_R23_upper
;;       
;;    endif
    
; ---------------------------------------------------------------------------    
; Pilyugin et al. (2001)
; ---------------------------------------------------------------------------    

    index = where((iabund.Z_R3 gt -900.0) and (iabund.Z_R23 gt -900.0),nindex)
    if (nindex ne 0L) then begin

; compute some ratios       
       
       R23 = 10^iabund[index].Z_R23
       R23_err = R23*alog(10.0)*iabund[index].Z_R23_err

       R3 = 10^iabund[index].Z_R3
       R3_err = R3*alog(10.0)*iabund[index].Z_R3_err
       
       P = R3/R23
       P_err = im_compute_error(R3,R3_err,R23,R23_err,/quotient)

; lower branch

       iabund[index].Z_12OH_P01_lower = 6.35 + 3.19*alog10(R23) - 1.74*alog10(R3)
       iabund[index].Z_12OH_P01_lower_err = sqrt( (3.19*R23_err/R23/alog(10))^2 + (1.74*R3_err/R3/alog(10))^2 )

; upper branch

       top = R23 + 54.2 + 59.45*P + 7.31*P^2
       top_err = sqrt( R23_err^2 + (59.45*P_err)^2 + (2.0*7.31*P*P_err)^2 )

       but = 6.07 + 6.71*P + 0.37*P^2 + 0.243*R23
       but_err = sqrt( (6.71*P_err)^2  + (2.0*0.37*P*P_err)^2 + (0.243*R23_err)^2 )

       P01 = top/but
       P01_err = im_compute_error(top,top_err,but,but_err,/quotient)
       
       iabund[index].Z_12OH_P01_upper = P01
       iabund[index].Z_12OH_P01_upper_err = P01_err

; compute the average of the two calibrations, to be used below

       iabund[index].Z_12OH_P01_average = (iabund[index].Z_12OH_P01_lower + iabund[index].Z_12OH_P01_upper) / 2.0
       iabund[index].Z_12OH_P01_average_err = sqrt((iabund[index].Z_12OH_P01_lower_err^2 + iabund[index].Z_12OH_P01_upper_err^2) / 2.0 )

; now, constrain the lower-branch abundances to be less than 7.95

       toobig = where(iabund[index].Z_12OH_P01_lower gt 7.95,ntoobig)
       if (ntoobig ne 0L) then begin
          iabund[index[toobig]].Z_12OH_P01_lower = -999.0
          iabund[index[toobig]].Z_12OH_P01_lower_err = -999.0
       endif
       
; and constrain the upper-branch abundances to be greater than 8.2

       toosmall = where(iabund[index].Z_12OH_P01_upper lt 8.2,ntoosmall)
       if (ntoosmall ne 0L) then begin
          iabund[index[toosmall]].Z_12OH_P01_upper = -999.0
          iabund[index[toosmall]].Z_12OH_P01_upper_err = -999.0
       endif
       
    endif
       
; ---------------------------------------------------------------------------
; use the procedure outlined in Melbourne & Salzer (2002) to figure
; out which Pilyugin branch we're on but use the new upper-branch
; calibration by Pilyugin (2002) instead of the Edmunds & Pagel
; calibration used by Melbourne & Salzer
; ---------------------------------------------------------------------------

    index = where((iabund.Z_nii_ha gt -900.0) and (iabund.Z_12OH_P01_upper gt -900.0),nindex)
    if (nindex ne 0L) then begin

; [N II/Ha < -1.3 --> lower branch
       
       lower = where((iabund[index].Z_nii_ha lt -1.3),nlower)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_P01_melbourne = iabund[index[lower]].Z_12OH_P01_lower
          iabund[index[lower]].Z_12OH_P01_melbourne_err = iabund[index[lower]].Z_12OH_P01_lower_err
       endif
       
; [N II/Ha > -1.0 --> upper branch
       
       upper = where((iabund[index].Z_nii_ha gt -1.0),nupper)
       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_P01_melbourne = iabund[index[upper]].Z_12OH_P01_upper
          iabund[index[upper]].Z_12OH_P01_melbourne_err = iabund[index[upper]].Z_12OH_P01_upper_err
       endif
       
    endif

; ---------------------------------------------------------------------------
; use O3N2 and N2 [Pettini] above and below 12+log(O/H)=8.5, respectively
; ---------------------------------------------------------------------------

    index = where((iabund.Z_12oh_n2_pettini gt -900.0) or (iabund.Z_12oh_O3N2_pettini gt -900.0),nindex)
    if (nindex ne 0L) then begin

       lower = where((iabund[index].Z_12oh_n2_pettini lt 8.5) and (iabund[index].Z_12oh_n2_pettini gt -900),nlower)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_O3N2_N2 = iabund[index[lower]].Z_12oh_n2_pettini
          iabund[index[lower]].Z_12OH_O3N2_N2_err = iabund[index[lower]].Z_12oh_n2_pettini_err
       endif
       
       upper = where((iabund[index].Z_12oh_n2_pettini ge 8.5) and (iabund[index].Z_12oh_n2_pettini gt -900) and $
         (iabund[index].Z_12oh_o3n2_pettini gt -900),nupper)
       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_O3N2_N2 = iabund[index[upper]].Z_12oh_O3N2_pettini
          iabund[index[upper]].Z_12OH_O3N2_N2_err = iabund[index[upper]].Z_12oh_O3N2_pettini_err
       endif
       
    endif

; ---------------------------------------------------------------------------
; use N2 [Pettini] to select the appropriate P01 calibration; we could
; use the average of the two calibrations between 12+log(O/H)=7.95 and
; 8.2 (and the commented code below allows you to do that), but this
; does not work well; so ignore the region between 7.95 and 8.2
; ---------------------------------------------------------------------------

    index = where((iabund.Z_12oh_n2_pettini gt -900.0) and ((iabund.Z_12oh_P01_upper gt -900.0) or $
      (iabund.Z_12oh_P01_lower gt -900.0)),nindex)
    if (nindex ne 0L) then begin

       lower = where((iabund[index].Z_12oh_n2_pettini lt 7.95) and (iabund[index].Z_12oh_P01_lower gt -900),nlower)
       if (nlower ne 0L) then begin
          iabund[index[lower]].Z_12OH_P01_N2 = iabund[index[lower]].Z_12oh_P01_lower
          iabund[index[lower]].Z_12OH_P01_N2_err = iabund[index[lower]].Z_12oh_P01_lower_err
       endif
       
       upper = where((iabund[index].Z_12oh_n2_pettini gt 8.2) and (iabund[index].Z_12oh_P01_upper gt -900),nupper)
       if (nupper ne 0L) then begin
          iabund[index[upper]].Z_12OH_P01_N2 = iabund[index[upper]].Z_12oh_P01_upper
          iabund[index[upper]].Z_12OH_P01_N2_err = iabund[index[upper]].Z_12oh_P01_upper_err
       endif
       
;      turn = where((iabund[index].Z_12oh_n2_pettini ge 7.95) and (iabund[index].Z_12oh_n2_pettini le 8.2) and $
;        (iabund[index].Z_12oh_P01_average gt -900),nturn)
;      if(nturn ne 0L) then begin
;         iabund[index[turn]].Z_12OH_P01_N2 = iabund[index[turn]].Z_12oh_P01_average
;         iabund[index[turn]].Z_12OH_P01_N2_err = iabund[index[turn]].Z_12oh_P01_average_err
;      endif

    endif

;;;; ---------------------------------------------------------------------------
;;;; use P01 on the upper and lower branches but the Melbourne & Salzer
;;;; (2002) calibration of [N II]/Ha in the turn-around region (see van
;;;; Zee et al. 1998)
;;;; ---------------------------------------------------------------------------
;;;
;;;    index = where((iabund.Z_12oh_n2 gt -900.0) or (iabund.Z_12oh_P gt -900.0),nindex)
;;;;   index = where((iabund.Z_12oh_n2 gt -900.0) or (iabund.Z_12oh_P gt -900.0) and $
;;;;     (iabund.Z_12OH_P01_upper gt -900.0),nindex)
;;;    if (nindex ne 0L) then begin
;;;
;;;;;; if both upper & lower branch P01 calibrations give the same branch
;;;;;; then use them!
;;;;;
;;;;;       upper = where((iabund[index].Z_12oh_p01_lower ge 8.2) and (iabund[index].Z_12OH_P01_upper ge 8.2),nupper)
;;;;;       if (nupper ne 0L) then begin
;;;;;          iabund[index[upper]].Z_12oh_P01 = iabund[index[upper]].Z_12oh_P01_upper
;;;;;          iabund[index[upper]].Z_12oh_P01_err = iabund[index[upper]].Z_12oh_P01_upper_err
;;;;;       endif
;;;;;       
;;;;;       lower = where((iabund[index].Z_12oh_p01_lower le 7.95) and (iabund[index].Z_12oh_p01_upper le 7.95),nlower)
;;;;;       if (nlower ne 0L) then begin
;;;;;          iabund[index[lower]].Z_12oh_P01 = iabund[index[lower]].Z_12oh_P01_lower
;;;;;          iabund[index[lower]].Z_12oh_P01_err = iabund[index[lower]].Z_12oh_P01_lower_err
;;;;;       endif
;;;
;;;; if we have N2 and P then use N2 otherwise use P
;;;
;;;       n2go = where((iabund[index].Z_12oh_n2 gt -900.0),nn2go)
;;;;      n2go = where((iabund[index].Z_12oh_n2 gt -900.0) and (iabund[index].Z_12oh_P01 lt -900.0),nn2go)
;;;       if (nn2go ne 0L) then begin
;;;
;;;; upper branch          
;;;          
;;;          upper = where((iabund[index[n2go]].Z_12oh_n2 ge 8.2) and $
;;;            (iabund[index[n2go]].Z_12OH_P01_upper gt -900.0),nupper)
;;;;         upper = where(iabund[index[n2go]].Z_12oh_n2 ge 8.2,nupper)
;;;
;;;          if (nupper ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nupper,format='(I0)')+$
;;;               ' upper-branch objects using N2-OH.'
;;;             iabund[index[n2go[upper]]].Z_12oh_P01 = iabund[index[n2go[upper]]].Z_12oh_P01_upper
;;;             iabund[index[n2go[upper]]].Z_12oh_P01_err = iabund[index[n2go[upper]]].Z_12oh_P01_upper_err
;;;          endif
;;;          
;;;; lower branch          
;;;
;;;          lower = where((iabund[index[n2go]].Z_12oh_n2 le 7.95) and $
;;;            (iabund[index[n2go]].Z_12OH_P01_lower gt -900.0),nlower)
;;;;         lower = where(iabund[index[n2go]].Z_12oh_n2 le 7.95,nlower)
;;;
;;;          if (nlower ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nlower,format='(I0)')+$
;;;               ' lower-branch objects using N2-OH.'
;;;             iabund[index[n2go[lower]]].Z_12oh_P01 = iabund[index[n2go[lower]]].Z_12oh_P01_lower
;;;             iabund[index[n2go[lower]]].Z_12oh_P01_err = iabund[index[n2go[lower]]].Z_12oh_P01_lower_err
;;;          endif
;;;
;;;; ambigious --> use N2
;;;
;;;          ambig = where((iabund[index[n2go]].Z_12oh_n2 gt 7.95) and (iabund[index[n2go]].Z_12oh_n2 lt 8.2),nambig)
;;;
;;;          if (nambig ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Applying N2-OH on '+string(nambig,format='(I0)')+' ambigious objects.'
;;;             iabund[index[n2go[ambig]]].Z_12oh_P01 = iabund[index[n2go[ambig]]].Z_12oh_n2
;;;             iabund[index[n2go[ambig]]].Z_12oh_P01_err = iabund[index[n2go[ambig]]].Z_12oh_n2_err
;;;          endif
;;;
;;;;;; ambigious --> use the lower/upper average
;;;;;
;;;;;          ambig = where((iabund[index[n2go]].Z_12oh_n2 gt 7.95) and (iabund[index[n2go]].Z_12oh_n2 lt 8.2),nambig)
;;;;;          if (nambig ne 0L) then begin
;;;;;             iabund[index[n2go[ambig]]].Z_12oh_P01 = $
;;;;;               (iabund[index[n2go[ambig]]].Z_12oh_P01_lower + iabund[index[n2go[ambig]]].Z_12oh_P01_upper) / 2.0
;;;;;             iabund[index[n2go[ambig]]].Z_12oh_P01_err = sqrt( iabund[index[n2go[ambig]]].Z_12oh_P01_lower_err^2 + $
;;;;;               iabund[index[n2go[ambig]]].Z_12oh_P01_upper_err^2)
;;;;;          endif
;;;
;;;       endif
;;;
;;;       Pgo = where((iabund[index].Z_12oh_P gt -900.0) and (iabund[index].Z_12oh_P01 lt -900.0),nPgo)
;;;       if (nPgo ne 0L) then begin
;;;
;;;          upper = where((iabund[index[Pgo]].Z_12oh_P ge 8.2) and $
;;;            (iabund[index[Pgo]].Z_12OH_P01_upper gt -900.0),nupper)
;;;
;;;          if (nupper ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nupper,format='(I0)')+$
;;;               ' upper-branch objects using P-OH.'
;;;             iabund[index[Pgo[upper]]].Z_12oh_P01 = iabund[index[Pgo[upper]]].Z_12oh_P01_upper
;;;             iabund[index[Pgo[upper]]].Z_12oh_P01_err = iabund[index[Pgo[upper]]].Z_12oh_P01_upper_err
;;;          endif
;;;          
;;;          lower = where((iabund[index[Pgo]].Z_12oh_P le 7.95) and $
;;;            (iabund[index[Pgo]].Z_12OH_P01_lower gt -900.0),nlower)
;;;
;;;          if (nlower ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Identified '+string(nlower,format='(I0)')+$
;;;               ' lower-branch objects using P-OH.'
;;;             iabund[index[Pgo[lower]]].Z_12oh_P01 = iabund[index[Pgo[lower]]].Z_12oh_P01_lower
;;;             iabund[index[Pgo[lower]]].Z_12oh_P01_err = iabund[index[Pgo[lower]]].Z_12oh_P01_lower_err
;;;          endif
;;;
;;;; ambigious --> use P
;;;
;;;          ambig = where((iabund[index[Pgo]].Z_12oh_P gt 7.95) and (iabund[index[Pgo]].Z_12oh_P lt 8.2),nambig)
;;;          if (nambig ne 0L) then begin
;;;             if (not keyword_set(silent)) then splog, 'Applying P-OH on '+string(nambig,format='(I0)')+' ambigious objects.'
;;;             iabund[index[Pgo[ambig]]].Z_12oh_P01 = iabund[index[Pgo[ambig]]].Z_12oh_P
;;;             iabund[index[Pgo[ambig]]].Z_12oh_P01_err = iabund[index[Pgo[ambig]]].Z_12oh_P_err
;;;          endif
;;;
;;;;;; ambigious --> use the lower/upper average
;;;;;
;;;;;          ambig = where((iabund[index[Pgo]].Z_12oh_P gt 7.95) and (iabund[index[Pgo]].Z_12oh_P lt 8.2),nambig)
;;;;;          if (nambig ne 0L) then begin
;;;;;             iabund[index[Pgo[ambig]]].Z_12oh_P01 = $
;;;;;               (iabund[index[Pgo[ambig]]].Z_12oh_P01_lower + iabund[index[Pgo[ambig]]].Z_12oh_P01_upper) / 2.0
;;;;;             iabund[index[Pgo[ambig]]].Z_12oh_P01_err = sqrt( iabund[index[Pgo[ambig]]].Z_12oh_P01_lower_err^2 + $
;;;;;               iabund[index[Pgo[ambig]]].Z_12oh_P01_upper_err^2)
;;;;;          endif
;;;
;;;       endif
;;;       
;;;;      upper = where(iabund[index].Z_nii_oii gt -0.7,nupper)
;;;;      if (nupper ne 0L) then iabund[index[upper]].Z_12OH_P01 = iabund[index[upper]].Z_12OH_P01_upper
;;;;      
;;;;      lower = where(iabund[index].Z_nii_oii lt -1.3,nlower)
;;;;      if (nlower ne 0L) then iabund[index[lower]].Z_12OH_P01 = iabund[index[lower]].Z_12OH_P01_lower
;;;;
;;;;      ambig = where(iabund[index].Z_12OH_P01 lt -900.0,nambig)
;;;;      if (nambig ne 0L) then begin
;;;;      
;;;;         upper_diff = abs(iabund[index[ambig]].Z_12oh_n2 - iabund[index[ambig]].Z_12OH_P01_upper)
;;;;         lower_diff = abs(iabund[index[ambig]].Z_12oh_n2 - iabund[index[ambig]].Z_12OH_P01_lower)
;;;;
;;;;         toobig = where((upper_diff ge 0.05) and (lower_diff ge 0.05),ntoobig,comp=okay,ncomp=nokay)
;;;;         if (ntoobig ne 0L) then iabund[index[ambig[toobig]]].Z_12OH_P01 = $
;;;;           iabund[index[ambig[toobig]]].Z_12oh_n2
;;;;
;;;;         if (nokay ne 0L) then begin
;;;;         
;;;;            upper = where((upper_diff[okay] lt lower_diff[okay]),nupper)
;;;;            if (nupper ne 0L) then iabund[index[ambig[okay[upper]]]].Z_12OH_P01 = $
;;;;              iabund[index[ambig[okay[upper]]]].Z_12OH_P01_upper
;;;;         
;;;;            lower = where((lower_diff[okay] lt upper_diff[okay]),nlower)
;;;;            if (nlower ne 0L) then iabund[index[ambig[okay[lower]]]].Z_12OH_P01 = $
;;;;              iabund[index[ambig[okay[lower]]]].Z_12OH_P01_lower
;;;;
;;;;         endif
;;;;
;;;;      endif
;;;
;;;    endif 

    abund[good_ebv] = iabund

return, abund
end
