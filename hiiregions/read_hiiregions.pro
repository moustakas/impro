;+
; NAME:
;   READ_HIIREGIONS()
;
; PURPOSE:
;   Read the HII-region database.  See BUILD_HIIREGIONS_DATABASE.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   nolog     - read the HII region structure where the fluxes
;               have not been "logarithm-ed"
;   hiiregion - only read in individual HII region data
;   hiigalaxy - only read in individual HII galaxy data
;   nosdss    - exclude data based on SDSS measurements
;   silent    - suppress messages to STDOUT
;   samplerefs - only return a small, representative sample of
;                references (MCCALL85, ZKH94, VANZEE98) 
;
; OUTPUTS:
;   hii - 
;
; OPTIONAL OUTPUTS:
;   hiined  - 
;   linefit - 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 July 05, U of A, written
;   jm04jul21uofa - added NOSDSS keyword
;   jm04dec17uofa - modified the directory structure
;   jm05may11uofa - treat individual HII regions separately from
;                   the HII galaxies; added HIIREGION and
;                   HIIGALAXY keywords; added HIINED optional
;                   output 
;   jm05may19uofa - added LINEFIT optional output
;   jm05aug19uofa - added SAMPLEREFS and NOKISS keywords
;   jm05sep27uofa - added LIMITEDREFS keyword
;   jm07dec19nyu  - added version number 
;
; Copyright (C) 2004-2005, John Moustakas
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

function read_hiiregions, hiined=hiined, linefit=linefit, bptplots=bptplots, $
  nolog=nolog, hiiregion=hiiregion, hiigalaxy=hiigalaxy, nosdss=nosdss, $
  silent=silent, samplerefs=samplerefs, limitedrefs=limitedrefs, nokiss=nokiss

    if (n_elements(bptplots) ne 0L) then nolog = 0L

    version = hiiregions_version()

    inpath = getenv('CATALOGS_DIR')+'/hiiregions/'
    if keyword_set(nolog) then $
      infile = 'hiiregions_database_nolog_'+version+'.fits.gz' else $
      infile = 'hiiregions_database_'+version+'.fits.gz'

    if (n_elements(silent) eq 0L) then splog, 'Reading '+inpath+infile
    hii = mrdfits(inpath+infile,1,/silent)
    if arg_present(linefit) then linefit = mrdfits(inpath+infile,2,/silent)

    if arg_present(hiined) then hiined = mrdfits(inpath+'hii_region_ned_basic_'+version+'.fits.gz',1,/silent)

    if keyword_set(hiiregion) then begin
       keep = where(hii.hiiregion,nkeep)
       if (nkeep ne 0L) then hii = hii[keep] else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif
    
    if keyword_set(hiigalaxy) then begin
       keep = where(hii.hiigalaxy,nkeep)
       if (nkeep ne 0L) then hii = hii[keep] else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif

    if keyword_set(nosdss) then begin
       keep = where((strmatch(hii.reference,'*Kniazev et al. 2004*',/fold) eq 0B) and $
         (strmatch(hii.reference,'*Izotov et al. 2006*',/fold) eq 0B),nkeep)
       if (nkeep ne 0L) then hii = hii[keep] else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif

    if keyword_set(nokiss) then begin
       keep = where(strmatch(hii.reference,'*Melbourne et al. 2004*',/fold) eq 0B,nkeep)
       if (nkeep ne 0L) then hii = hii[keep] else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif

    if keyword_set(samplerefs) then begin
       keep = where($
         (strmatch(hii.reference,'*McCall et al. 1985*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Izotov, Thuan, & Lipovetsky 1997*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Izotov, Thuan, & Lipovetsky 1994*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Izotov & Thuan 1998*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Zaritsky et al. 1994*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*van Zee et al. 1998*',/fold) eq 1B),nkeep)
       if (nkeep ne 0L) then begin
          hii = hii[keep]
          if arg_present(linefit) then linefit = linefit[keep]
       endif else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif
    
    if keyword_set(limitedrefs) then begin
       keep = where($
         (strmatch(hii.reference,'*McCall et al. 1985*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Izotov & Thuan 1998*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*Zaritsky et al. 1994*',/fold) eq 1B) or $
         (strmatch(hii.reference,'*van Zee et al. 1998*',/fold) eq 1B),nkeep)
       if (nkeep ne 0L) then hii = hii[keep] else begin
          print, 'No matching data.'
          return, hii
       endelse
    endif
    
; Fix this later!    

;   splog, 'Warning: Excluding Lee et al.!!'
;   keep = where(strmatch(hii.reference,'*Lee et al. 2003*',/fold) eq 0B,nkeep)
;   if (nkeep ne 0L) then hii = hii[keep] else begin
;      print, 'No matching data.'
;      return, hii
;   endelse

    if keyword_set(bptplots) then begin

       models = kewley_bpt_lines()

       im_window, 0, xratio=0.5, /square

; ---------------------------------------------------------------------------       

       good = where((hii.oiii_5007_h_beta gt -900.0) and (hii.nii_6584_h_alpha gt -900.0))
       Te = where(hii[good].ZT_log12oh gt -900.0,comp=nonTe)

       xTe = hii[good[Te]].nii_6584_h_alpha
       xerrTe = hii[good[Te]].nii_6584_h_alpha_err
       yTe = hii[good[Te]].oiii_5007_h_beta
       yerrTe = hii[good[Te]].oiii_5007_h_beta_err

       if (nonTe ne 0L) then begin
          xnonTe = hii[good[nonTe]].nii_6584_h_alpha
          xerrnonTe = hii[good[nonTe]].nii_6584_h_alpha_err
          ynonTe = hii[good[nonTe]].oiii_5007_h_beta
          yerrnonTe = hii[good[nonTe]].oiii_5007_h_beta_err
       endif

       Z = hii[good[Te]].ZT_log12oh
       
       plotsym, 0, 0.5
       djs_plot, xnonTe, ynonTe, ps=8, xsty=3, ysty=3, xrange=[-2.7,0.5], $
         yrange=[-2.0,1.5], xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
         color='red', xtitle='log ([N II] / H\alpha)', ytitle='log ([O III] / H\beta)'

       plotsym, 8, 0.5
       djs_oplot, xTe, yTe, ps=8, color='grey'

;      ploterror, x, y, xerr, yerr, ps=3, xsty=3, ysty=3, xrange=[-2.6,1.0], $
;    yrange=[-2.0,1.5], xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0

       oplot, models.x_nii, models.y_nii, line=0, thick=2.0
       cc = get_kbrd(1)

; ---------------------------------------------------------------------------       
       
       good = where((hii.oiii_5007_h_beta gt -900.0) and (hii.sii_h_alpha gt -900.0))
       Te = where(hii[good].ZT_log12oh gt -900.0,comp=nonTe)

       xTe = hii[good[Te]].sii_h_alpha
       xerrTe = hii[good[Te]].sii_h_alpha_err

       xnonTe = hii[good[nonTe]].sii_h_alpha
       xerrnonTe = hii[good[nonTe]].sii_h_alpha_err

       yTe = hii[good[Te]].oiii_5007_h_beta
       yerrTe = hii[good[Te]].oiii_5007_h_beta_err

       ynonTe = hii[good[nonTe]].oiii_5007_h_beta
       yerrnonTe = hii[good[nonTe]].oiii_5007_h_beta_err

       Z = hii[good[Te]].ZT_log12oh
       
       plotsym, 0, 0.5
       djs_plot, xnonTe, ynonTe, ps=8, xsty=3, ysty=3, xrange=[-2.4,0.5], $
         yrange=[-2.0,1.5], xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
         color='red', xtitle='log ([S II] / H\alpha)', ytitle='log ([O III] / H\beta)'

       plotsym, 8, 0.5
       djs_oplot, xTe, yTe, ps=8, color='grey'

       oplot, models.x_sii, models.y_sii, line=0, thick=2.0
       cc = get_kbrd(1)
       
; ---------------------------------------------------------------------------       

    endif
    
return, hii
end
