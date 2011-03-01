; ###########################################################################
; THIS VERSION ALSO CLASSIFIES BASED ON [SII] and [OI]!
; ###########################################################################
;+
; NAME:
;   ICLASSIFICATION_EXPANDED()
;
; PURPOSE:
;   Classify galaxies into starburst or AGN based on their
;   emission-line ratios.
;
; INPUTS:
;   linenodust - ispec-style data structure
;
; OPTIONAL INPUTS:
;   snrcut_class - minimum S/N to apply to each emission line (default
;                  1.0) 
;   snrcut_nii   - optionally adopt a different S/N cut for [NII]
;                  (default 1.0)
;   chi2cut      - maximum chi2 to apply to each emission line
;                  (default 1E8, i.e., none)
;   sigmacut     - maximum sigma cut to apply to each emission line
;                  (default 500 km/s)
;   niihacut     - classify objects with log([NII]/Ha)>NIIHACUT as AGN
;                  (default -0.4)
;
; KEYWORD PARAMETERS:
;   kauffmann - use the Kauffmann et al. (2003) 
;   doplot    - generate a plot of each classification diagram and
;               wait for a key-stroke
;   silent    - suppress messages to STDOUT
;
; OUTPUTS:
;   bptclass - output data structure with the classification results
;
; OPTIONAL OUTPUTS:
;   ratios - emission-line ratios used in the classification
;
; COMMENTS:
;   We use observed line fluxes to do the classification.  Note
;   that this routine does not treat upper limits.
;
; TODO:
;   Identify galaxies in the transition region.  Also identify
;   galaxies that are classified based on fewer than three
;   diagrams.  
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 June 13, U of A
;   jm04may13uofa - added new physical criteria for classifying
;                   ambigious objects
;   jm05jun09uofa - require a 1-sigma detection on each line-ratio 
;   jm05nov02uofa - error checking for GALAXY tag in LINE 
;   jm06feb08uofa - compute the MIXTURE class (see below); on
;                   error, return BPTCLASS rather than -1
;   jm06oct14nyu  - added SNRCUT_CLASS and NIIHACUT optional
;                   inputs 
;   jm07aug10nyu  - minor code improvements
;   jm08feb10nyu  - return the errors in the RATIOS structure;
;                   added SNRCUT_NII optional input; improved
;                   plotting 
;   jm08feb25nyu  - added BPT_D and BPT_PHI output
;   jm08aug26nyu  - added additional quality cuts, requiring that the
;                   flux *error* be positive; documentation tweaks;
;                   added CHI2CUT and SIGMACUT optional inputs
;
; Copyright (C) 2003-2008, John Moustakas
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

function iclassification_expanded, line, ratios=ratios, kauffmann=kauffmann, $
  snrcut_class=snrcut_class, snrcut_nii=snrcut_nii, chi2cut=chi2cut, $
  sigmacut=sigmacut, niihacut=niihacut, doplot=doplot, silent=silent

    nspec = n_elements(line)
    if nspec eq 0L then begin
       doc_library, 'iclassification'
       return, -1L
    endif

    if (n_elements(snrcut_class) eq 0L) then snrcut_class = 1.0
    if (n_elements(snrcut_nii) eq 0L) then snrcut_nii = snrcut_class
    if (n_elements(chi2cut) eq 0L) then chi2cut = 5.0
    if (n_elements(sigmacut) eq 0L) then sigmacut = 500.0
    if (n_elements(niihacut) eq 0L) then niihacut = -0.4
    
    if (not keyword_set(silent)) then $
      splog, 'S/N > '+string(snrcut_class,format='(G0.0)')+'; '+$
      'S/N [N II] > '+string(snrcut_nii,format='(G0.0)')+'; '+$
      'Chi2 < '+string(chi2cut,format='(G0.0)')+'; '+$
      'Sigma < '+string(sigmacut,format='(G0.0)')+' km/s'

; initialize the output structure    
    
    bptclass = {$
      bpt_id:                      0L, $
      bpt_galaxy:                 '?', $
      bpt_nii_kauffmann_class:     0L, $ ; 0 = no classification; 2 = AGN; 4 = HII
      bpt_nii_kewley_class:        0L, $
      bpt_nii_class:               0L, $
      bpt_sii_class:               0L, $
      bpt_oi_class:                0L, $
      bpt_d:                   -999.0, $ ; Kauffmann "D" value (distance from starburst sequence)
      bpt_phi:                 -999.0, $ ; Kauffmann "Phi" (degrees, measured positive clockwise)
      bpt_class:                  '?', $ ; AGN, HII, AMB (ambiguous), or NO DATA
      bpt_class_quality:           0L, $ ; classification based on CLASS_QUALITY diagrams (0, 1, 2, or 3)
      bpt_pure_nii_class:         '?', $ ; classification based on the BPT diagram and [N II]/H-alpha
      bpt_nii_mixture_class:      '?', $ ; classification based purely on [N II]/H-alpha -->
      bpt_ebv:                  0.0}     ; --> but based on a mix of the Kewley and Kauffmann curves 
    bptclass = replicate(bptclass,nspec)

    bptclass.bpt_id = lindgen(nspec)
    if tag_exist(line,'GALAXY') then bptclass.bpt_galaxy = line.galaxy

; construct the line ratio structure necessary to classify each object 

    ratios = {$
      galaxy:          '', $
      nii_ha:      -999.0, $
      nii_ha_err:  -999.0, $
      sii_ha:      -999.0, $
      sii_ha_err:  -999.0, $
      oi_ha:       -999.0, $
      oi_ha_err:   -999.0, $
      oiii_hb:     -999.0, $
      oiii_hb_err: -999.0, $
      class:           ''}
    ratios = replicate(ratios,nspec)
    if tag_exist(line,'GALAXY') then ratios.galaxy = line.galaxy
    
; verify that LINE has sufficient information to do the classification
; and fill RATIOS

    tags = tag_names(line)

    if tag_exist(line,'H_ALPHA') then begin
       if tag_exist(line,'H_ALPHA_CHI2') then hachi2 = line.h_alpha_chi2 else hachi2 = fltarr(nspec)
       if tag_exist(line,'H_ALPHA_SIGMA') then hasigma = line.h_alpha_sigma[0] else hasigma = fltarr(nspec)
       hagood = where((line.h_alpha[0]/line.h_alpha[1] ge snrcut_class) and $
         (line.h_alpha[1] gt 0.0) and (hachi2 lt chi2cut) and (hasigma lt sigmacut),nhagood) 
       if nhagood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'Insufficient H-alpha line fluxes'
          return, bptclass
       endif
    endif else begin
       if (not keyword_set(silent)) then splog, 'BPT classification requires H-alpha line fluxes'
       return, bptclass
    endelse

    if tag_exist(line,'H_BETA') then begin
       if tag_exist(line,'H_BETA_CHI2') then hbchi2 = line.h_beta_chi2 else hbchi2 = fltarr(nspec)
       if tag_exist(line,'H_BETA_SIGMA') then hbsigma = line.h_beta_sigma[0] else hbsigma = fltarr(nspec)
       hbgood = where((line.h_beta[0]/line.h_beta[1] ge snrcut_class) and $
         (line.h_beta[1] gt 0.0) and (hbchi2 lt chi2cut) and (hbsigma lt sigmacut),nhbgood) 
       if nhbgood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'Insufficient H-beta line fluxes'
          return, bptclass
       endif
    endif else begin
       if (not keyword_set(silent)) then splog, 'BPT classification requires H-beta line fluxes'
       return, bptclass
    endelse

    if tag_exist(line,'OIII_5007') then begin
       if tag_exist(line,'OIII_5007_CHI2') then oiiichi2 = line.oiii_5007_chi2 else oiiichi2 = fltarr(nspec)
       if tag_exist(line,'OIII_5007_SIGMA') then oiiisigma = line.oiii_5007_sigma[0] else oiiisigma = fltarr(nspec)
       oiiigood = where((line.oiii_5007[0]/line.oiii_5007[1] ge snrcut_class) and $
         (line.oiii_5007[1] gt 0.0) and (oiiichi2 lt chi2cut) and (oiiisigma lt sigmacut),noiiigood) 
       if noiiigood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'Insufficient [O III] 5007 line fluxes'
          return, bptclass
       endif
    endif else begin
       if (not keyword_set(silent)) then splog, 'BPT classification requires [O III] 5007 line fluxes'
       return, bptclass
    endelse

    if (nhbgood ne 0L) and (noiiigood ne 0L) then begin
       oiiihbgood = where((line.oiii_5007[0]/line.oiii_5007[1] ge snrcut_class) and $
         (line.oiii_5007[1] gt 0.0) and (oiiichi2 lt chi2cut) and (oiiisigma lt sigmacut) and $
         (line.h_beta[0]/line.h_beta[1] ge snrcut_class) and (line.h_beta[1] gt 0.0) and $
         (hbchi2 lt chi2cut) and (hbsigma lt sigmacut),noiiihbgood)

       if (noiiihbgood eq 0L) then begin
          if (not keyword_set(silent)) then splog, 'Insufficient data to form [O III]/H-beta line ratio'
          return, bptclass       
       endif else begin
          lineratio, line[oiiihbgood], 'OIII_5007', 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=0.0
          ratios[oiiihbgood].oiii_hb = x
          ratios[oiiihbgood].oiii_hb_err = xerr
;         ratios[oiiihbgood].oiii_hb = alog10(line[oiiihbgood].oiii_5007[0]/line[oiiihbgood].h_beta[0]) ; [O III]/H-beta
       endelse
    endif
          
; read the theoretical classification lines    
    
    models_kauffmann = kewley_bpt_lines(/kauffmann)
    models_kewley = kewley_bpt_lines()

    if keyword_set(kauffmann) then models = models_kauffmann else models = models_kewley

; ---------------------------------------------------------------------------
; NII/Ha    
; ---------------------------------------------------------------------------

    niipossible = 1
    if tag_exist(line,'NII_6584') then begin
       if tag_exist(line,'NII_6584_CHI2') then niichi2 = line.nii_6584_chi2 else niichi2 = fltarr(nspec)
       if tag_exist(line,'NII_6584_SIGMA') then niisigma = line.nii_6584_sigma[0] else niisigma = fltarr(nspec)

       niigood = where((line.nii_6584[0]/line.nii_6584[1] ge snrcut_class) and $
         (line.nii_6584[1] gt 0.0) and (niichi2 lt chi2cut) and (niisigma lt sigmacut),nniigood)
       if nniigood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'No well-measured [N II] line fluxes'
          niipossible = 0
       endif else begin
          if (nhagood ne 0L) and (nniigood ne 0L) then begin
             niihagood = where((line.nii_6584[0]/line.nii_6584[1] ge snrcut_nii) and $ ; NOTE!
               (line.nii_6584[1] gt 0.0) and (niichi2 lt chi2cut) and (niisigma lt sigmacut) and $
               (line.h_alpha[0]/line.h_alpha[1] ge snrcut_class) and (line.h_alpha[1] gt 0.0) and $
               (hachi2 lt chi2cut) and (hasigma lt sigmacut),nniihagood)

             if (nniihagood ne 0L) then begin
                lineratio, line[niihagood], 'NII_6584', 'H_ALPHA', '', '', $
                  x, xerr, index=index, nindex=nindex, snrcut=0.0
                ratios[niihagood].nii_ha = x
                ratios[niihagood].nii_ha_err = xerr
;               ratios[niihagood].nii_ha = alog10(line[niihagood].nii_6584[0]/line[niihagood].h_alpha[0]) ; [N II]/H-alpha
             endif else begin
                if (not keyword_set(silent)) then splog, 'No well-measured [N II]/H-alpha line ratios'
                niipossible = 0
             endelse 
          endif else niipossible = 0
             
       endelse 

       if niipossible then begin ; do the classification

          good = cmset_op(niihagood,'AND',oiiihbgood)
          if good[0] ne -1L then begin

             bptclass[good].bpt_d = bpt_dphi(ratios[good].nii_ha,ratios[good].oiii_hb,phi=bpt_phi)
             bptclass[good].bpt_phi = bpt_phi
             
; -------------------------
; Kauffmann classification             
; -------------------------
             
             data = kewley_bpt_lines(ratio_nii=ratios[good].nii_ha,/nii,/kauffmann)
             agn_indx = ratios[good].oiii_hb gt data.y_nii

             w = where(ratios[good].nii_ha gt max(models.x_nii),nw)
             if nw ne 0L then agn_indx[w] = 1

             agn = where(agn_indx,nagn,comp=hii,ncomp=nhii)
             if nagn ne 0L then bptclass[good[agn]].bpt_nii_kauffmann_class = 2
             if nhii ne 0L then bptclass[good[hii]].bpt_nii_kauffmann_class = 4

; -------------------------
; Kewley classification             
; -------------------------
             
             data = kewley_bpt_lines(ratio_nii=ratios[good].nii_ha,/nii)
             agn_indx = ratios[good].oiii_hb gt data.y_nii

             w = where(ratios[good].nii_ha gt max(models.x_nii),nw)
             if nw ne 0L then agn_indx[w] = 1

             agn = where(agn_indx,nagn,comp=hii,ncomp=nhii)
             if nagn ne 0L then bptclass[good[agn]].bpt_nii_kewley_class = 2
             if nhii ne 0L then bptclass[good[hii]].bpt_nii_kewley_class = 4

; -------------------------

          endif 

       endif  

    endif else begin

       if (not keyword_set(silent)) then splog, 'No [N II]/H-alpha classification possible'
       niipossible = 0

    endelse

    if keyword_set(kauffmann) then $
      bptclass.bpt_nii_class = bptclass.bpt_nii_kauffmann_class else $
      bptclass.bpt_nii_class = bptclass.bpt_nii_kewley_class
    
; ---------------------------------------------------------------------------
; SII/Ha
; ---------------------------------------------------------------------------

    siipossible = 1
    if tag_exist(line,'SII_6716') and tag_exist(line,'SII_6731') then begin
       if tag_exist(line,'SII_6716_CHI2') then sii6716chi2 = line.sii_6716_chi2 else sii6716chi2 = fltarr(nspec)
       if tag_exist(line,'SII_6731_CHI2') then sii6731chi2 = line.sii_6731_chi2 else sii6731chi2 = fltarr(nspec)
       if tag_exist(line,'SII_6716_SIGMA') then sii6716sigma = line.sii_6716_sigma[0] else sii6716sigma = fltarr(nspec)
       if tag_exist(line,'SII_6731_SIGMA') then sii6731sigma = line.sii_6731_sigma[0] else sii6731sigma = fltarr(nspec)
      
       siigood = where((line.sii_6716[0]/line.sii_6716[1] ge snrcut_class) and $
         (line.sii_6716[1] gt 0.0) and (sii6716chi2 lt chi2cut) and (sii6716sigma lt sigmacut) and $
         (line.sii_6731[0]/line.sii_6731[1] ge snrcut_class) and $
         (line.sii_6731[1] gt 0.0) and (sii6731chi2 lt chi2cut) and (sii6731sigma lt sigmacut),nsiigood)
       if nsiigood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'No well-measured [S II] line fluxes'
          siipossible = 0
       endif else begin

          if (nhagood ne 0L) and (nsiigood ne 0L) then begin
             siihagood = where((line.sii_6716[0]/line.sii_6716[1] ge snrcut_class) and $
               (line.sii_6716[1] gt 0.0) and (sii6716chi2 lt chi2cut) and (sii6716sigma lt sigmacut) and $
               (line.sii_6731[0]/line.sii_6731[1] ge snrcut_class) and $
               (line.sii_6731[1] gt 0.0) and (sii6731chi2 lt chi2cut) and (sii6731sigma lt sigmacut) and $
               (line.h_alpha[0]/line.h_alpha[1] ge snrcut_class) and (line.h_alpha[1] gt 0.0) and $
               (hachi2 lt chi2cut) and (hasigma lt sigmacut),nsiihagood)

             if (nsiihagood ne 0L) then begin
                lineratio, line[siihagood], ['SII_6716','SII_6731'], 'H_ALPHA', '', '', $
                  x, xerr, index=index, nindex=nindex, snrcut=0.0
                ratios[siihagood].sii_ha = x
                ratios[siihagood].sii_ha_err = xerr
;               ratios[siihagood].sii_ha = alog10((line[siihagood].sii_6716[0]+$
;                 line[siihagood].sii_6731[0])/line[siihagood].h_alpha[0]) ; [S II]/H-alpha
             endif else begin
                if (not keyword_set(silent)) then splog, 'No well-measured [S II]/H-alpha line ratios'
                siipossible = 0
             endelse 
          endif else siipossible = 0
             
       endelse 

       if siipossible then begin ; do the classification

          good = cmset_op(siihagood,'AND',oiiihbgood)
          if good[0] ne -1L then begin
          
             data = kewley_bpt_lines(ratio_sii=ratios[good].sii_ha,/sii)
             agn_indx = ratios[good].oiii_hb gt data.y_sii

; this WHERE statement catches AGNs located beyond the well-defined
; region of the hyperbola
             
             w = where(ratios[good].sii_ha gt max(models.x_sii),nw)
             if nw ne 0L then agn_indx[w] = 1

             agn = where(agn_indx,nagn,comp=hii,ncomp=nhii)
             if nagn ne 0L then bptclass[good[agn]].bpt_sii_class = 2
             if nhii ne 0L then bptclass[good[hii]].bpt_sii_class = 4

          endif 

       endif

    endif else begin

       if (not keyword_set(silent)) then splog, 'No [S II]/H-alpha classification possible'
       siipossible = 0

    endelse

; ---------------------------------------------------------------------------
; OI/Ha
; ---------------------------------------------------------------------------

    oipossible = 1
    if tag_exist(line,'OI_6300') then begin
       if tag_exist(line,'OI_6300_CHI2') then oichi2 = line.oi_6300_chi2 else oichi2 = fltarr(nspec)
       if tag_exist(line,'OI_6300_SIGMA') then oisigma = line.oi_6300_sigma[0] else oisigma = fltarr(nspec)

       oigood = where((line.oi_6300[0]/line.oi_6300[1] ge snrcut_class) and $
         (line.oi_6300[1] gt 0.0) and (oichi2 lt chi2cut) and (oisigma lt sigmacut),noigood)
       if noigood eq 0L then begin
          if (not keyword_set(silent)) then splog, 'No well-measured [O I] line fluxes'
          oipossible = 0
       endif else begin

          if (nhagood ne 0L) and (noigood ne 0L) then begin
             oihagood = where((line.oi_6300[0]/line.oi_6300[1] ge snrcut_class) and $
               (line.oi_6300[1] gt 0.0) and (oichi2 lt chi2cut) and (oisigma lt sigmacut) and $
               (line.h_alpha[0]/line.h_alpha[1] ge snrcut_class) and (line.h_alpha[1] gt 0.0) and $
               (hachi2 lt chi2cut) and (hasigma lt sigmacut),noihagood)
             if (noihagood ne 0L) then begin
                lineratio, line[oihagood], 'OI_6300', 'H_ALPHA', '', '', $
                  x, xerr, index=index, nindex=nindex, snrcut=0.0
                ratios[oihagood].oi_ha = x
                ratios[oihagood].oi_ha_err = xerr
;               ratios[oihagood].oi_ha = alog10(line[oihagood].oi_6300[0]/line[oihagood].h_alpha[0]) ; [O I]/H-alpha
             endif else begin
                if (not keyword_set(silent)) then splog, 'No well-measured [O I]/H-alpha line ratios'
                oipossible = 0
             endelse 
          endif else oipossible = 0
             
       endelse 

       if oipossible then begin ; do the classification

          good = cmset_op(oihagood,'AND',oiiihbgood)
          if good[0] ne -1L then begin
          
             data = kewley_bpt_lines(ratio_oi=ratios[good].oi_ha,/oi)
             agn_indx = ratios[good].oiii_hb gt data.y_oi

; this WHERE statement catches AGNs located beyond the well-defined
; region of the hyperbola
             
             w = where(ratios[good].oi_ha gt max(models.x_oi),nw)
             if nw ne 0L then agn_indx[w] = 1

             agn = where(agn_indx,nagn,comp=hii,ncomp=nhii)
             if nagn ne 0L then bptclass[good[agn]].bpt_oi_class = 2
             if nhii ne 0L then bptclass[good[hii]].bpt_oi_class = 4

          endif 

       endif

    endif else begin

       if (not keyword_set(silent)) then splog, 'No [O I]/H-alpha classification possible'
       oipossible = 0

    endelse

; do the final classification based on all three diagrams

    for k = 0L, nspec-1L do begin

       bpt = [bptclass[k].bpt_nii_class,bptclass[k].bpt_sii_class,bptclass[k].bpt_oi_class]
       noclass = where(bpt eq 0.0,nnoclass,comp=class,ncomp=nclass)
       bptclass[k].bpt_class_quality = nclass

       type = 'Unknown'
       if nclass ne 0L then begin
          
          case total(bpt) of
             nclass*2.0: type = 'AGN'
             nclass*4.0: type = 'HII'
             else: type = 'AMB'
          endcase

       endif

       bptclass[k].bpt_class = type
       
    endfor

; add a classification based purely on the [N II]/H-alpha ratio;
; objects with no classification have poorly measured [O III] or
; H-beta fluxes

    bptclass.bpt_pure_nii_class = 'Unknown'
    
    hii = where(bptclass.bpt_nii_class eq 4L,nhii)
    if (nhii ne 0L) then bptclass[hii].bpt_pure_nii_class = 'HII'
    
    agn = where(bptclass.bpt_nii_class eq 2L,nagn)
    if (nagn ne 0L) then bptclass[agn].bpt_pure_nii_class = 'AGN'

    unknown = where((bptclass.bpt_nii_class eq 0L) and (ratios.nii_ha gt -999.0),nunknown)
    if (nunknown ne 0L) then begin
;      niceprint, line[unknown].oiii_5007[1], line[unknown].nii_6584[1], $
;        line[unknown].h_beta[1], line[unknown].h_alpha[1]
       hii = where(ratios[unknown].nii_ha lt niihacut,nhii,comp=agn,ncomp=nagn)
       if (nhii ne 0L) then bptclass[unknown[hii]].bpt_pure_nii_class = 'HII'
       if (nagn ne 0L) then bptclass[unknown[agn]].bpt_pure_nii_class = 'AGN'
    endif

; determine the mixture class    

    bptclass.bpt_nii_mixture_class = 'Unknown'
    
    hii = where((bptclass.bpt_nii_kauffmann_class eq 4L) and (bptclass.bpt_nii_kewley_class eq 4L),nhii)
    if (nhii ne 0L) then bptclass[hii].bpt_nii_mixture_class = 'HII'
    
    agn = where((bptclass.bpt_nii_kauffmann_class eq 2L) and (bptclass.bpt_nii_kewley_class eq 2L),nagn)
    if (nagn ne 0L) then bptclass[agn].bpt_nii_mixture_class = 'AGN'
    
    mix = where(((bptclass.bpt_nii_kauffmann_class eq 2L) and (bptclass.bpt_nii_kewley_class eq 4L)) or $
      ((bptclass.bpt_nii_kauffmann_class eq 2L) and (bptclass.bpt_nii_kewley_class eq 4L)),nmix)
    if (nmix ne 0L) then bptclass[mix].bpt_nii_mixture_class = 'HII/AGN'

; ---------------------------------------------------------------------------    
; attempt to classify some of the ambigious galaxies using additional
; constraints:
; ---------------------------------------------------------------------------    

; because the [O I]/Ha ratio is much more sensitive to shocks in the
; integrated spectra of galaxies, if an object is an AGN in [O I]/Ha
; but HII in the [N II]/Ha and [S II]/Ha diagrams then re-classify it
; as an HII galaxy

    amb = where(bptclass.bpt_class eq 'AMB',namb)
    if namb ne 0L then begin

       reclass = where((bptclass[amb].bpt_nii_class eq 4L) and (bptclass[amb].bpt_sii_class eq 4L) and $
         (bptclass[amb].bpt_oi_class eq 2L),nreclass)
       if (nreclass ne 0L) then begin
          if (not keyword_set(silent)) then splog, 'Relegating the [O I]/Ha class in '+string(nreclass,format='(I0)')+'/'+$
            string(namb,format='(I0)')+' galaxies: AMB --> HII'
          bptclass[amb[reclass]].bpt_class = 'HII'
          bptclass[amb[reclass]].bpt_class_quality = bptclass[amb[reclass]].bpt_class_quality - 1L
       endif

    endif

; if there are additional ambigious objects classify them based on the
; [N II]/Ha ratio because of the shock-sensitivity of the [O I]/Ha
; ratio and the density-sensitivity of [S II]/Ha

    amb = where(bptclass.bpt_class eq 'AMB',namb)
    if namb ne 0L then begin

       hii = where(bptclass[amb].bpt_nii_class eq 4L,nhii,comp=agn,ncomp=nagn)
       if (nhii ne 0L) then begin
          if (not keyword_set(silent)) then splog, 'Relegating the [S II]/Ha, [O I]/Ha classes in '+string(nhii,format='(I0)')+'/'+$
            string(namb,format='(I0)')+' galaxies: AMB --> HII'
          bptclass[amb[hii]].bpt_class = 'HII'
          bptclass[amb[hii]].bpt_class_quality = 1L
       endif
       if (nagn ne 0L) then begin
          if (not keyword_set(silent)) then splog, 'Relegating the [S II]/Ha, [O I]/Ha classes in '+string(nagn,format='(I0)')+'/'+$n
            string(namb,format='(I0)')+' galaxies: AMB --> AGN'
          bptclass[amb[agn]].bpt_class = 'AGN'
          bptclass[amb[agn]].bpt_class_quality = 1L
       endif

    endif

; finally attempt to classify "Unknown objects" based entirely on a
; measurement of the [N II]/Ha ratio; most of these objects should
; have poorly measured [O III]/Hb ratios, which may bias us against
; high-metallicity objects (likely AGN)
    
    unknown = where(bptclass.bpt_class eq 'Unknown',nunknown)
    if nunknown ne 0L then begin

       good = where(ratios[unknown].nii_ha gt -900.0,ngood)

       if ngood ne 0L then begin
          hii = where(ratios[unknown[good]].nii_ha lt niihacut,nhii,comp=agn,ncomp=nagn)
          if nhii ne 0L then begin
             if (not keyword_set(silent)) then splog, 'Reclassifying based on [N II]/Ha '+string(nhii,format='(I0)')+'/'+$
               string(nunknown,format='(I0)')+' galaxies: Unknown --> HII'
             bptclass[unknown[good[hii]]].bpt_class = 'HII'
             bptclass[unknown[good[hii]]].bpt_class_quality = 1L
          endif
          if nagn ne 0L then begin
             if (not keyword_set(silent)) then splog, 'Reclassifying based on [N II]/Ha '+string(nagn,format='(I0)')+'/'+$
               string(nunknown,format='(I0)')+' galaxies: Unknown --> AGN'
             bptclass[unknown[good[agn]]].bpt_class = 'AGN'
             bptclass[unknown[good[agn]]].bpt_class_quality = 1L
          endif
       endif 
    
    endif 

    ratios.class = bptclass.bpt_class
    
;  struct_print, bptclass

    agn = where(bptclass.bpt_class eq 'AGN',nagn)
    hii = where(bptclass.bpt_class eq 'HII',nhii)
    amb = where(bptclass.bpt_class eq 'AMB',namb)
    insuf = where(bptclass.bpt_class eq 'Unknown',ninsuf)
    
    if (not keyword_set(silent)) then splog, format='("Found ",I0,"/",I0," HII galaxies, ",I0,"/",I0," AGN, "'+$
      ',I0,"/",I0," ambigious classifications, and ",I0,"/",I0,'+$
           '" unknown (insufficient data).")', $
      nhii, nspec, nagn, nspec, namb, nspec, ninsuf, nspec

    agn = where(bptclass.bpt_pure_nii_class eq 'AGN',nagn)
    hii = where(bptclass.bpt_pure_nii_class eq 'HII',nhii)
    insuf = where(bptclass.bpt_pure_nii_class eq 'Unknown',ninsuf)
    
    if (not keyword_set(silent)) then splog, format='("Found ",I0,"/",I0," HII galaxies, ",I0,"/",I0," AGN, and "'+$
      ',I0,"/",I0," unknown (insufficient data) galaxies based on [N II]/Ha only.")', $
      nhii, nspec, nagn, nspec, ninsuf, nspec

; ---------------------------------------------------------------------------    
; generate the three BPT diagrams
; ---------------------------------------------------------------------------    

    if keyword_set(doplot) then begin

;      window, 0, xs=550, ys=550

       if niipossible then begin
          
          agn = where(bptclass.bpt_nii_class eq 2,nagn)
          hii = where(bptclass.bpt_nii_class eq 4,nhii)
          insuf = where(bptclass.bpt_nii_class eq 0,ninsuf)
          
          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, thick=2.0, $
            xrange=[-2.0,1.0], yrange=[-1.4,1.5], xthick=2.0, ythick=2.0, $
            charsize=1.8, charthick=2.0, xtitle='log ([N II] \lambda6584/H\alpha)', $
            ytitle='log ([O III] \lambda5007/H\beta)'

          if nagn ne 0L then begin
;            oploterror, ratios[agn].nii_ha, ratios[agn].oiii_hb, $
;              ratios[agn].nii_ha_err, ratios[agn].oiii_hb_err, ps=5, $
;              color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
             if nagn eq 1L then plots, ratios[agn].nii_ha, ratios[agn].oiii_hb, ps=5 else $
               djs_oplot, ratios[agn].nii_ha, ratios[agn].oiii_hb, ps=6, sym=0.25, color='red'
          endif 
          if nhii ne 0L then begin
;            oploterror, ratios[hii].nii_ha, ratios[hii].oiii_hb, $
;              ratios[hii].nii_ha_err, ratios[hii].oiii_hb_err, ps=4, $
;              color=djs_icolor('blue'), errcolor=djs_icolor('blue')
             if nhii eq 1L then plots, ratios[hii].nii_ha, ratios[hii].oiii_hb, ps=2 else $
               djs_oplot, ratios[hii].nii_ha, ratios[hii].oiii_hb, ps=6, sym=0.25, color='blue'
          endif 
;;          if namb ne 0L then begin
;;             oploterror, ratios[amb].nii_ha, ratios[amb].oiii_hb, $
;;               ratios[amb].nii_ha_err, ratios[amb].oiii_hb_err, ps=4, $
;;               color=djs_icolor('red'), errcolor=djs_icolor('red')
;;;            if namb eq 1L then plots, ratios[amb].nii_ha, ratios[amb].oiii_hb, ps=2 else $
;;;              djs_oplot, ratios[amb].nii_ha, ratios[amb].oiii_hb, ps=4, color='green'
;;          endif 

          oplot, models_kewley.x_nii, models_kewley.y_nii, line=2, thick=2.0
;         oplot, models_kewley.x_nii, models_kewley.y_nii_upper, line=2, thick=2.0
;         oplot, models_kewley.x_nii, models_kewley.y_nii_lower, line=2, thick=2.0
          oplot, models_kauffmann.x_nii, models_kauffmann.y_nii, line=0, thick=2.0

          splog, '[N II]/Ha: press any key to continue'
          cc = get_kbrd(1)

       endif

;;       if siipossible then begin
;;          
;;          djs_plot, models.x_sii, models.y_sii, xsty=3, ysty=3, thick=2.0, $
;;            xrange=[-1.5,0.5], yrange=[-1.4,1.5], xthick=2.0, ythick=2.0, $
;;            charsize=1.8, charthick=2.0, xtitle='log ([S II]/H\alpha)', $
;;            ytitle='log ([O III]/H\beta)'
;;          oplot, models.x_sii, models.y_sii_upper, line=2, thick=2.0
;;          oplot, models.x_sii, models.y_sii_lower, line=2, thick=2.0
;;          
;;          if nagn ne 0L then begin
;;             if nagn eq 1L then plots, ratios[agn].sii_ha, ratios[agn].oiii_hb, ps=5 else $
;;               djs_oplot, ratios[agn].sii_ha, ratios[agn].oiii_hb, ps=4, color='red'
;;          endif
;;          if nhii ne 0L then begin
;;             if nhii eq 1L then plots, ratios[hii].sii_ha, ratios[hii].oiii_hb, ps=2 else $
;;               djs_oplot, ratios[hii].sii_ha, ratios[hii].oiii_hb, ps=4, color='blue'
;;          endif
;;          if namb ne 0L then begin
;;             if namb eq 1L then plots, ratios[amb].sii_ha, ratios[amb].oiii_hb, ps=2 else $
;;               djs_oplot, ratios[amb].sii_ha, ratios[amb].oiii_hb, ps=4, color='green'
;;          endif
;;          splog, '[S II]/Ha: press any key to continue.'
;;          cc = get_kbrd(1)       
;;
;;       endif
;;
;;       if oipossible then begin
;;          
;;          djs_plot, models.x_oi, models.y_oi, xsty=3, ysty=3, thick=2.0, $
;;            xrange=[-2.5,0.0], yrange=[-1.4,1.5], xthick=2.0, ythick=2.0, $
;;            charsize=1.8, charthick=2.0, xtitle='log ([O I]/H\alpha)', $
;;            ytitle='log ([O III]/H\beta)'
;;          oplot, models.x_oi, models.y_oi_upper, line=2, thick=2.0
;;          oplot, models.x_oi, models.y_oi_lower, line=2, thick=2.0
;;          
;;          if nagn ne 0L then begin
;;             if nagn eq 1L then plots, ratios[agn].oi_ha, ratios[agn].oiii_hb, ps=5 else $
;;               djs_oplot, ratios[agn].oi_ha, ratios[agn].oiii_hb, ps=4, color='red'
;;          endif
;;          if nhii ne 0L then begin
;;             if nhii eq 1L then plots, ratios[hii].oi_ha, ratios[hii].oiii_hb, ps=2 else $
;;               djs_oplot, ratios[hii].oi_ha, ratios[hii].oiii_hb, ps=4, color='blue'
;;          endif
;;          if namb ne 0L then begin
;;             if namb eq 1L then plots, ratios[amb].oi_ha, ratios[amb].oiii_hb, ps=2 else $
;;               djs_oplot, ratios[amb].oi_ha, ratios[amb].oiii_hb, ps=4, color='green'
;;          endif
;;          splog, '[O I]/Ha : press any key to continue.'
;;          cc = get_kbrd(1)
;;
;;       endif
          
    endif 
    
return, bptclass
end
