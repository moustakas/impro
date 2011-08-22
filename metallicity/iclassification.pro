;+
; NAME:
;   ICLASSIFICATION()
;
; PURPOSE:
;   Classify galaxies into starburst or AGN based on their
;   emission-line ratios.
;
; INPUTS:
;   line - ispec-style emission-line data structure
;
; OPTIONAL INPUTS:
;   snrcut_class - minimum S/N to apply to each emission line; can be
;     a scalar or a 4-element array corresponding to each emission
;     line: H-alpha, H-beta, [OIII] 5007, and [NII] 6584, in that
;     order (default is 1.0 for all lines)
;   niihacut - classify objects with log([NII]/Ha)>NIIHACUT as AGN
;     (default -0.3) 
;
; KEYWORD PARAMETERS:
;   noniihaclass - do not classify using *just* the [NII]/Ha ratio
;     (NIIHACUT is ignored)
;   doplot - generate a diagnostic plot
;   silent - suppress messages to STDOUT
;
; OUTPUTS:
;   bptclass - bit mask indicating the cumulative classification
;     results (see COMMENTS)
;
; OPTIONAL OUTPUTS:
;   ratios - emission-line ratios used in the classification as well
;     as the final classifications (NOCLASS, UNKNOWN, AGN, SF, and
;     SF/AGN) 
;
; COMMENTS:
;   Write this!
;
; EXAMPLES:
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
;   jm08sep26nyu  - bit of a major rewrite; streamlined to just use
;     the [NII]/Ha vs [OIII]/Hb diagram; now uses bits to store the
;     classifications 
;   jm10jan21ucsd - removed CHI2CUT and SIGMACUT: the user is
;     responsible for inputting reasonable line-measurements; also
;     removed SNRCUT_NII and allowed SNRCUT_CLASS to be a 4-element
;     array for each of the four emission lines 
;
; Copyright (C) 2003-2009, John Moustakas
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

function iclassification, line, ratios=ratios, snrcut_class=snrcut_class, $
  niihacut=niihacut, noniihaclass=noniihaclass, doplot=doplot, silent=silent

    nspec = n_elements(line)
    if nspec eq 0L then begin
       doc_library, 'iclassification'
       return, -1L
    endif

; parse SNRCUT_CLASS    
    if (n_elements(snrcut_class) eq 0) then begin
       snrcut_ha = 1.0
       snrcut_hb = 1.0
       snrcut_oiii = 1.0
       snrcut_nii = 1.0
    endif else begin
       nsnrcut = n_elements(snrcut_class)
       if (nsnrcut ne 1) and (nsnrcut ne 4) then begin
          splog, 'SNRCUT_CLASS must be a 1- or 4-element array'
          return, -1
       endif
       if (nsnrcut eq 1) then begin
          if (keyword_set(silent) eq 0) then splog, 'S/N > '+$
            string(snrcut_class,format='(G0.0)')
          snrcut_ha = snrcut_class
          snrcut_hb = snrcut_class
          snrcut_oiii = snrcut_class
          snrcut_nii = snrcut_class
       endif else begin
          snrcut_ha = snrcut_class[0]
          snrcut_hb = snrcut_class[1]
          snrcut_oiii = snrcut_class[2]
          snrcut_nii = snrcut_class[3]
          if (keyword_set(silent) eq 0) then begin
             splog, 'S/N H-alpha > '+string(snrcut_ha,format='(G0.0)')
             splog, 'S/N H-beta  > '+string(snrcut_hb,format='(G0.0)')
             splog, 'S/N [O III] > '+string(snrcut_oiii,format='(G0.0)')
             splog, 'S/N [N II]  > '+string(snrcut_nii,format='(G0.0)')
          endif
       endelse
    endelse

    if (n_elements(niihacut) eq 0) then niihacut = -0.3
    
; initialize the output bit array and the structure of line-ratios
; necessary to classify each object  
    bptclass = lonarr(nspec)

    ratios = {$
      galaxy:          '', $
      bpt_d:       -999.0, $ ; Kauffmann "D" value (distance from starburst sequence)
      bpt_phi:     -999.0, $ ; Kauffmann "Phi" (degrees, measured positive clockwise)
      nii_ha:      -999.0, $
      nii_ha_err:  -999.0, $
;     sii_ha:      -999.0, $
;     sii_ha_err:  -999.0, $
;     oi_ha:       -999.0, $
;     oi_ha_err:   -999.0, $
      oiii_hb:     -999.0, $
      oiii_hb_err: -999.0, $
      final_class:     ''}
    ratios = replicate(ratios,nspec)
    if tag_exist(line,'GALAXY') then ratios.galaxy = line.galaxy
    
; verify that LINE has sufficient information to do the classification
; and fill RATIOS
    tags = tag_names(line)

    if tag_exist(line,'H_ALPHA') then begin
       hagood = where((line.h_alpha[0]/line.h_alpha[1] ge snrcut_ha) and $
         (line.h_alpha[1] gt 0.0),nhagood) 
;      if nhagood eq 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Insufficient H-alpha line fluxes'
;         bptclass = bptclass or class_flagval('CLASS','NOCLASS')
;         return, bptclass
;      endif
    endif else begin
       if (keyword_set(silent) eq 0) then splog, 'BPT classification requires H-alpha line fluxes'
       bptclass = bptclass or class_flagval('CLASS','NOCLASS')
       return, bptclass
    endelse

    if tag_exist(line,'H_BETA') then begin
       hbgood = where((line.h_beta[0]/line.h_beta[1] ge snrcut_hb) and $
         (line.h_beta[1] gt 0.0),nhbgood) 
;      if nhbgood eq 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Insufficient H-beta line fluxes'
;         bptclass = bptclass or class_flagval('CLASS','NOCLASS')
;         return, bptclass
;      endif
    endif else begin
       if (keyword_set(silent) eq 0) then splog, 'BPT classification requires H-beta line fluxes'
       bptclass = bptclass or class_flagval('CLASS','NOCLASS')
       return, bptclass
    endelse

    if tag_exist(line,'OIII_5007') then begin
       oiiigood = where((line.oiii_5007[0]/line.oiii_5007[1] ge snrcut_oiii) and $
         (line.oiii_5007[1] gt 0.0),noiiigood) 
;      if noiiigood eq 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Insufficient [O III] 5007 line fluxes'
;         bptclass = bptclass or class_flagval('CLASS','NOCLASS')
;         return, bptclass
;      endif
    endif else begin
       if (keyword_set(silent) eq 0) then splog, 'BPT classification requires [O III] 5007 line fluxes'
       bptclass = bptclass or class_flagval('CLASS','NOCLASS')
       return, bptclass
    endelse

    if tag_exist(line,'NII_6584') then begin
       niigood = where((line.nii_6584[0]/line.nii_6584[1] ge snrcut_nii) and $
         (line.nii_6584[1] gt 0.0),nniigood) 
;      if nniigood eq 0L then begin
;         if (keyword_set(silent) eq 0) then splog, 'Insufficient [N II] 6584 line fluxes'
;         bptclass = bptclass or class_flagval('CLASS','NOCLASS')
;         return, bptclass
;      endif
    endif else begin
       if (keyword_set(silent) eq 0) then splog, 'BPT classification requires [N II] 6584 line fluxes'
       bptclass = bptclass or class_flagval('CLASS','NOCLASS')
       return, bptclass
    endelse

; ---------------------------------------------------------------------------
; [OIII]/Hb
    if (nhbgood ne 0L) and (noiiigood ne 0L) then begin
       oiiihbgood = where((line.oiii_5007[0]/line.oiii_5007[1] ge snrcut_oiii) and $
         (line.oiii_5007[1] gt 0.0) and (line.h_beta[0]/line.h_beta[1] ge snrcut_hb) and $
         (line.h_beta[1] gt 0.0),noiiihbgood)
       if (noiiihbgood ne 0L) then begin
          lineratio, line[oiiihbgood], 'OIII_5007', 'H_BETA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=0.0
          ratios[oiiihbgood].oiii_hb = x
          ratios[oiiihbgood].oiii_hb_err = xerr
;         ratios[oiiihbgood].oiii_hb = alog10(line[oiiihbgood].oiii_5007[0]/line[oiiihbgood].h_beta[0]) ; [O III]/H-beta
       endif
    endif else noiiihbgood = 0L
          
; ---------------------------------------------------------------------------
; [NII]/Ha    
    if (nhagood ne 0L) and (nniigood ne 0L) then begin
       niihagood = where((line.nii_6584[0]/line.nii_6584[1] ge snrcut_nii) and $
         (line.nii_6584[1] gt 0.0) and (line.h_alpha[0]/line.h_alpha[1] ge snrcut_ha) and $
         (line.h_alpha[1] gt 0.0),nniihagood)
       if (nniihagood ne 0L) then begin
          lineratio, line[niihagood], 'NII_6584', 'H_ALPHA', '', '', $
            x, xerr, index=index, nindex=nindex, snrcut=0.0
          ratios[niihagood].nii_ha = x
          ratios[niihagood].nii_ha_err = xerr
       endif
    endif else nniihagood = 0L

; read the classification curves and then classify!
    models_kauffmann = kewley_bpt_lines(/kauffmann)
    models_kewley = kewley_bpt_lines()

    if (noiiihbgood gt 0L) and (nniihagood gt 0L) then begin
       good = cmset_op(niihagood,'AND',oiiihbgood)
       if (good[0] ne -1L) then begin

          ratios[good].bpt_d = bpt_dphi(ratios[good].nii_ha,ratios[good].oiii_hb,phi=bpt_phi)
          ratios[good].bpt_phi = bpt_phi
             
; -------------------------
; Kauffmann classification             
          data = kewley_bpt_lines(ratio_nii=ratios[good].nii_ha,/nii,/kauffmann)
          agn_indx = ratios[good].oiii_hb gt data.y_nii

          w = where(ratios[good].nii_ha gt max(models_kauffmann.x_nii),nw)
          if nw ne 0L then agn_indx[w] = 1

          agn = where(agn_indx,nagn,comp=sf,ncomp=nsf)
          if nagn ne 0L then bptclass[good[agn]] = bptclass[good[agn]] or $
            class_flagval('CLASS','KAUFFMANN_AGN')
          if nsf ne 0L then bptclass[good[sf]] = bptclass[good[sf]] or $
            class_flagval('CLASS','KAUFFMANN_SF')

; -------------------------
; Kewley classification             
          data = kewley_bpt_lines(ratio_nii=ratios[good].nii_ha,/nii)
          agn_indx = ratios[good].oiii_hb gt data.y_nii

          w = where(ratios[good].nii_ha gt max(models_kewley.x_nii),nw)
          if nw ne 0L then agn_indx[w] = 1

          agn = where(agn_indx,nagn,comp=sf,ncomp=nsf)
          if nagn ne 0L then bptclass[good[agn]] = bptclass[good[agn]] or $
            class_flagval('CLASS','KEWLEY_AGN')
          if nsf ne 0L then bptclass[good[sf]] = bptclass[good[sf]] or $
            class_flagval('CLASS','KEWLEY_SF')
       endif 
    endif

 ; add a classification based purely on the [N II]/H-alpha ratio
    if (keyword_set(noniihaclass) eq 0) then begin
       niiagn = where((ratios.nii_ha gt -999.0) and (ratios.nii_ha gt niihacut),nniiagn)
       if (nniiagn ne 0L) then bptclass[niiagn] = bptclass[niiagn] or $
         class_flagval('CLASS','NII_AGN')

       niisf = where((ratios.nii_ha gt -999.0) and (ratios.nii_ha lt niihacut),nniisf)
       if (nniisf ne 0L) then bptclass[niisf] = bptclass[niisf] or $
         class_flagval('CLASS','NII_SF')
    endif

; assign an unknown classification to anything that hasn't yet
; been classified
    unk = where((bptclass eq 0),nunk)
    if (nunk ne 0L) then bptclass[unk] = bptclass[unk] or $
      class_flagval('CLASS','UNKNOWN')
    
; generate the final classifications: KEWLEY_AGN=AGN; KAUFFMANN_SF=SF;
; between KAUFFMANN_SF and KEWLEY_AGN=SF/AGN
    flag = where(bptclass eq 0) & if flag[0] ne -1 then message, 'This is bad'
    
    noclass = where((bptclass and class_flagval('CLASS','NOCLASS')),nnoclass)
    unk = where((bptclass and class_flagval('CLASS','UNKNOWN')),nunk)
    agn = where((bptclass and class_flagval('CLASS','KEWLEY_AGN')),nagn)
    sf = where((bptclass and class_flagval('CLASS','KAUFFMANN_SF')),nsf)
    sfagn = where(((bptclass and class_flagval('CLASS','KAUFFMANN_AGN')) ne 0) and $
      ((bptclass and class_flagval('CLASS','KEWLEY_SF')) ne 0),nsfagn)
;   help, bptclass, nagn, nsf, nsfagn, nnoclass, nunk

; what about the NII_AGN/NII_SF classes?; for now give them separate
; classifications, but only for things that have not been previously
; classified (note 'EQ')
    niiagn = where((bptclass eq class_flagval('CLASS','NII_AGN')),nniiagn)
    niisf = where((bptclass eq class_flagval('CLASS','NII_SF')),nniisf)

    if (nnoclass ne 0L) then ratios[noclass].final_class = 'NOCLASS'
    if (nunk ne 0L) then ratios[unk].final_class = 'UNKNOWN'

; store the final classes
    ncanclass = 0
    if (nagn ne 0L) then begin
       ratios[agn].final_class = 'AGN'
       ncanclass = ncanclass + nagn
    endif
    if (nsf ne 0L) then begin
       ratios[sf].final_class = 'SF'
       ncanclass = ncanclass + nsf
    endif
    if (nsfagn ne 0L) then begin
       ratios[sfagn].final_class = 'SF/AGN'
       ncanclass = ncanclass + nsfagn
    endif
    if (nniiagn ne 0L) then begin
       ratios[niiagn].final_class = 'NII_AGN'
       ncanclass = ncanclass + nniiagn
    endif
    if (nniisf ne 0L) then begin
       ratios[niisf].final_class = 'NII_SF'
       ncanclass = ncanclass + nniisf
    endif

    flag = where(strtrim(ratios.final_class,2) eq '')
    if flag[0] ne -1 then message, 'This is bad'
    
    if (keyword_set(silent) eq 0) then begin
       splog, 'Not classifiable: '+string(nnoclass,format='(I0)')+'/'+string(nspec,format='(I0)')+' ('+$
         strtrim(string(100.0*nnoclass/nspec,format='(F12.1)'),2)+'%)'
       splog, 'Unclassified: '+string(nunk,format='(I0)')+'/'+string(nspec,format='(I0)')+' ('+$
         strtrim(string(100.0*nunk/nspec,format='(F12.1)'),2)+'%)'
       splog, 'Classified  : '+string(ncanclass,format='(I0)')+'/'+string(nspec,format='(I0)')+' ('+$
         strtrim(string(100.0*ncanclass/nspec,format='(F12.1)'),2)+'%)'
       splog, ' BPT-AGN    : '+string(nagn,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
         strtrim(string(100.0*nagn/ncanclass,format='(F12.1)'),2)+'%)'
       splog, ' BPT-SF     : '+string(nsf,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
         strtrim(string(100.0*nsf/ncanclass,format='(F12.1)'),2)+'%)'
       splog, ' BPT-SF/AGN : '+string(nsfagn,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
         strtrim(string(100.0*nsfagn/ncanclass,format='(F12.1)'),2)+'%)'
       if (nniiagn ne 0L) then $
         splog, ' NII-AGN    : '+string(nniiagn,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
         strtrim(string(100.0*nniiagn/ncanclass,format='(F12.1)'),2)+'%)'
       if (nniisf ne 0L) then $
         splog, ' NII-SF     : '+string(nniisf,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
         strtrim(string(100.0*nniisf/ncanclass,format='(F12.1)'),2)+'%)'
    endif

; ---------------------------------------------------------------------------    
; generate a BPT diagram
    if keyword_set(doplot) then begin
       djs_plot, [0], [0], /nodata, xsty=3, ysty=3, thick=2.0, $
         xrange=[-2.0,1.0], yrange=[-1.4,1.5], xthick=2.0, ythick=2.0, $
         charsize=1.8, charthick=2.0, xtitle='log ([N II] \lambda6584/H\alpha)', $
         ytitle='log ([O III] \lambda5007/H\beta)'

       if nagn ne 0L then begin
;         oploterror, ratios[agn].nii_ha, ratios[agn].oiii_hb, $
;           ratios[agn].nii_ha_err, ratios[agn].oiii_hb_err, ps=5, $
;           color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
          if nagn eq 1L then plots, ratios[agn].nii_ha, ratios[agn].oiii_hb, ps=5 else $
            djs_oplot, ratios[agn].nii_ha, ratios[agn].oiii_hb, ps=6, sym=0.25, color='red'
       endif 
       if nsf ne 0L then begin
;         oploterror, ratios[sf].nii_ha, ratios[sf].oiii_hb, $
;           ratios[sf].nii_ha_err, ratios[sf].oiii_hb_err, ps=4, $
;           color=djs_icolor('blue'), errcolor=djs_icolor('blue')
          if nsf eq 1L then plots, ratios[sf].nii_ha, ratios[sf].oiii_hb, ps=2 else $
            djs_oplot, ratios[sf].nii_ha, ratios[sf].oiii_hb, ps=6, sym=0.25, color='blue'
       endif 
       if nsfagn ne 0L then begin
;         oploterror, ratios[sfagn].nii_ha, ratios[sfagn].oiii_hb, $
;           ratios[sfagn].nii_ha_err, ratios[sfagn].oiii_hb_err, ps=4, $
;           color=djs_icolor('blue'), errcolor=djs_icolor('blue')
          if nsfagn eq 1L then plots, ratios[sfagn].nii_ha, ratios[sfagn].oiii_hb, ps=2 else $
            djs_oplot, ratios[sfagn].nii_ha, ratios[sfagn].oiii_hb, ps=6, sym=0.25, color='dark green'
       endif 

       oplot, models_kewley.x_nii, models_kewley.y_nii, line=2, thick=2.0
;      oplot, models_kewley.x_nii, models_kewley.y_nii_upper, line=2, thick=2.0
;      oplot, models_kewley.x_nii, models_kewley.y_nii_lower, line=2, thick=2.0
       oplot, models_kauffmann.x_nii, models_kauffmann.y_nii, line=0, thick=2.0

;      splog, '[N II]/Ha: press any key to continue'
;      cc = get_kbrd(1)
    endif 
    
    if keyword_set(name) then begin
    endif
        
;;    agn = where(bptclass.bpt_pure_nii_class eq 'AGN',nagn)
;;    sf = where(bptclass.bpt_pure_nii_class eq 'SF',nsf)
;;    insuf = where(bptclass.bpt_pure_nii_class eq 'Unknown',ninsuf)
;;    
;;    if (not keyword_set(silent)) then splog, format='("Found ",I0,"/",I0," SF galaxies, ",I0,"/",I0," AGN, and "'+$
;;      ',I0,"/",I0," unknown (insufficient data) galaxies based on [N II]/Ha only.")', $
;;      nsf, nspec, nagn, nspec, ninsuf, nspec

return, bptclass
end
