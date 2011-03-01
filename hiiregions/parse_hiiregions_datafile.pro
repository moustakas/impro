;+
; NAME:
;   PARSE_HIIREGIONS_DATAFILE
;
; PURPOSE:
;   Parse one of the standardized HII-region datafiles (see
;   BUILD_HIIREGIONS_FLUXES).
;
; INPUTS: 
;   datafile - 
;   reference - 
;   texref - 
;
; KEYWORD PARAMETERS: 
;
;
; OUTPUTS: 
;   data - 
;
; OPTIONAL OUTPUTS:
;   data_nolog - 
;   linefit - 
;
; COMMENTS:
;   I assume that there are no good fluxes (flux ne -999) without good
;   flux errors. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, ??
;   jm10feb22ucsd - documented and updated
;
; Copyright (C) 2010, John Moustakas
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

function parse_hiiregions_datafile, datafile, reference, $
  texref, linefit=linefit, data_nolog=data_nolog

; branching ratios
    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    nratio = branch.n_ii  ; 3.054
    sratio = branch.s_iii ; 2.480

    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio
    scor = 1.0+1.0/sratio

    ndatafile = n_elements(datafile)
    for j = 0L, ndatafile-1L do begin

; read the data file       
       if (file_test(datafile[j],/regular) eq 0L) then begin
          splog, 'Data file '+datafile[j]+' not found!'
          stop
       endif
       
       hiidata1 = rsex(datafile[j])
       nobject = n_elements(hiidata1)

; some datafiles were typed in (e.g., 2003_kennicutt.sex) and others
; were automatically parsed (e.g., 2006_izotov.sex); this trick makes
; sure everything is living in the same structure; note that the older
; datafiles used TOIII instead of T4363 so copy that separately 
       hiidata = init_one_hiiregion_structure(nobject)
       struct_assign, hiidata1, hiidata, /nozero

       if (tag_exist(hiidata1,'t4363') eq 0B) and (tag_exist(hiidata1,'toiii') eq 1B) then begin
          hiidata.t4363     = hiidata1.toiii
          hiidata.t4363_err = hiidata1.toiii_err
       endif
       
; ---------------------------------------------------------------------------
; initialize the output data structure and store the miscellaneous HII
; region properties
       d = init_hiiregions_structure(nobject)

       d.hii_galaxy = hiidata.hii_galaxy
       d.hii_region = hiidata.hii_region

       d.hii_raoffset   = hiidata.raoffset
       d.hii_deoffset   = hiidata.deoffset
       d.hii_lit_radius = hiidata.radius
       d.hii_lit_rr25   = hiidata.rr25
       d.ewhb           = hiidata.ewhb
       d.fha            = hiidata.fha

       d.reference     = reference[j]
       d.texref        = texref[j]

; ---------------------------------------------------------------------------
; line-fluxes that exist for all data files (historically)       
       d.oii_h_beta           = hiidata.oii_3727
       d.oii_h_beta_err       = hiidata.oii_3727_err
       d.oiii_4363_h_beta     = hiidata.oiii_4363
       d.oiii_4363_h_beta_err = hiidata.oiii_4363_err
       d.oiii_5007_h_beta     = hiidata.oiii_5007
       d.oiii_5007_h_beta_err = hiidata.oiii_5007_err
       d.nii_6584_h_beta      = hiidata.nii_6584
       d.nii_6584_h_beta_err  = hiidata.nii_6584_err
       d.sii_6716_h_beta      = hiidata.sii_6716
       d.sii_6716_h_beta_err  = hiidata.sii_6716_err
       d.sii_6731_h_beta      = hiidata.sii_6731
       d.sii_6731_h_beta_err  = hiidata.sii_6731_err

; ---------------------------------------------------------------------------
; literature electron temperatures and abundances; note that we do not
; check for T[NII] as we *always* assume that T[NII]=T[OII]  
       d.lit_log12oh_te      = hiidata.log12oh
       d.lit_log12oh_te_err  = hiidata.log12oh_err

       g = where(hiidata.t4363 gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].lit_t4363     = 1E4*hiidata[g].t4363
          d[g].lit_t4363_err = 1E4*hiidata[g].t4363_err
       endif
       g = where(hiidata.t5755 gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].lit_t5755     = 1E4*hiidata[g].t5755
          d[g].lit_t5755_err = 1E4*hiidata[g].t5755_err
       endif
       g = where(hiidata.t6312 gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].lit_t6312     = 1E4*hiidata[g].t6312
          d[g].lit_t6312_err = 1E4*hiidata[g].t6312_err
       endif
       g = where(hiidata.t7325 gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].lit_t7325     = 1E4*hiidata[g].t7325
          d[g].lit_t7325_err = 1E4*hiidata[g].t7325_err
       endif

;; older datafiles use TOIII instead of T4363
;       if (tag_exist(hiidata,'t4363') eq 0B) and (tag_exist(hiidata,'toiii') eq 1B) then begin
;          g = where(hiidata.toiii gt -900.0,ng)
;          if (ng ne 0L) then begin
;             d[g].lit_t4363     = 1E4*hiidata[g].toiii
;             d[g].lit_t4363_err = 1E4*hiidata[g].toiii_err
;          endif
;       endif

       ha     = hiidata.ha
       ha_err = hiidata.ha_err

;      haflag = where((ha lt -900.0),nhaflag)
;      if (nhaflag ne 0L) then splog, 'No Ha for region: '+$
;        strjoin(hiidata[haflag].hii_region,',')+$
;        ' ('+file_basename(datafile[j])+')'
       
       noha = where((ha lt -900.0),nnoha)
       if (nnoha ne 0L) then begin
          te = where((d[noha].lit_t4363 gt -900.0),nte,comp=note,ncomp=nnote)
          if (nte ne 0L) then begin
             ha[noha[te]] = return_tbalmer(d[noha[te]].lit_t4363,/HaHb)
             ha_err[noha[te]] = 0.05*ha[noha[te]]
          endif
          if (nnote ne 0L) then begin
             ha[noha[note]] = return_tbalmer(1D4,/HaHb)
             ha_err[noha[note]] = 0.05*ha[noha[note]]
          endif
       endif

       haflag = where((ha lt -900.0),nhaflag)
       if (nhaflag ne 0L) then message, 'Fix me!'
       
       d.h_alpha_h_beta       = ha
       d.h_alpha_h_beta_err   = ha_err

       ha_good = ha gt -900.0

; ---------------------------------------------------------------------------
; [O II] 3727
       g = where((d.oii_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].oii_h_alpha     = d[g].oii_h_beta/ha[g]
          d[g].oii_h_alpha_err = d[g].oii_h_beta_err
;         d[g].oii_h_alpha_err = im_compute_error(d[g].oii_h_beta,d[g].oii_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

; ---------------------------------------------------------------------------
; [O II] 7325
       d.oii_7325_h_beta     = hiidata.oii_7325
       d.oii_7325_h_beta_err = hiidata.oii_7325_err
       g = where((d.oii_7325_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].oii_7325_h_alpha     = d[g].oii_7325_h_beta/ha[g]
          d[g].oii_7325_h_alpha_err = d[g].oii_7325_h_beta_err
;         d[g].oii_7325_h_alpha_err = im_compute_error(d[g].oii_7325_h_beta,$
;           d[g].oii_7325_h_beta_err,ha[g],ha_err[g],/quotient)
       endif
       
; ---------------------------------------------------------------------------
; [O III] 4363
       g = where((d.oiii_4363_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].oiii_4363_h_alpha = d[g].oiii_4363_h_beta/ha[g]
          d[g].oiii_4363_h_alpha_err = d[g].oiii_4363_h_beta_err
;         d[g].oiii_4363_h_alpha_err = im_compute_error(d[g].oiii_4363_h_beta,$
;           d[g].oiii_4363_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

; ---------------------------------------------------------------------------
; [O III] 4959, 5007
       g = where((d.oiii_5007_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].oiii_5007_h_alpha = d[g].oiii_5007_h_beta/ha[g]
          d[g].oiii_5007_h_alpha_err = d[g].oiii_5007_h_beta_err
;         d[g].oiii_5007_h_alpha_err = im_compute_error(d[g].oiii_5007_h_beta,$
;           d[g].oiii_5007_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

       g = where(d.oiii_5007_h_beta gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].oiii_4959_h_beta     = d[g].oiii_5007_h_beta/oratio
          d[g].oiii_4959_h_beta_err = d[g].oiii_5007_h_beta_err/oratio
       endif

       g = where(d.oiii_5007_h_alpha gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].oiii_4959_h_alpha     = d[g].oiii_5007_h_alpha/oratio
          d[g].oiii_4959_h_alpha_err = d[g].oiii_5007_h_alpha_err/oratio
       endif

       g = where(d.oiii_5007_h_beta gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].oiii_h_beta     = d[g].oiii_4959_h_beta+d[g].oiii_5007_h_beta
          d[g].oiii_h_beta_err = sqrt(d[g].oiii_4959_h_beta_err^2+d[g].oiii_5007_h_beta_err^2)
       endif

       g = where(d.oiii_5007_h_alpha gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].oiii_h_alpha     = d[g].oiii_4959_h_alpha+d[g].oiii_5007_h_alpha
          d[g].oiii_h_alpha_err = sqrt(d[g].oiii_4959_h_alpha_err^2+d[g].oiii_5007_h_alpha_err^2)
       endif

; ---------------------------------------------------------------------------
; [N II] 5755
       d.nii_5755_h_beta           = hiidata.nii_5755
       d.nii_5755_h_beta_err       = hiidata.nii_5755_err
       g = where((d.nii_5755_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].nii_5755_h_alpha     = d[g].nii_5755_h_beta/ha[g]
          d[g].nii_5755_h_alpha_err = d[g].nii_5755_h_beta_err
;         d[g].nii_5755_h_alpha_err = im_compute_error(d[g].nii_5755_h_beta,$
;           d[g].nii_5755_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

; ---------------------------------------------------------------------------
; [N II] 6548, 6584
       g = where((d.nii_6584_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].nii_6584_h_alpha = d[g].nii_6584_h_beta/ha[g]
          d[g].nii_6584_h_alpha_err = d[g].nii_6584_h_beta_err
;         d[g].nii_6584_h_alpha_err = im_compute_error(d[g].nii_6584_h_beta,$
;           d[g].nii_6584_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

       g = where(d.nii_6584_h_beta gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].nii_6548_h_beta     = d[g].nii_6584_h_beta/nratio
          d[g].nii_6548_h_beta_err = d[g].nii_6584_h_beta_err/nratio
       endif
       g = where(d.nii_6584_h_alpha gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].nii_6548_h_alpha     = d[g].nii_6584_h_alpha/nratio
          d[g].nii_6548_h_alpha_err = d[g].nii_6584_h_alpha_err/nratio
       endif

       g = where(d.nii_6584_h_beta gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].nii_h_beta     = d[g].nii_6548_h_beta+d[g].nii_6584_h_beta
          d[g].nii_h_beta_err = sqrt(d[g].nii_6548_h_beta_err^2+d[g].nii_6584_h_beta_err^2)
       endif
       g = where(d.nii_6584_h_alpha gt -900.0,ng)
       if (ng ne 0L) then begin
          d[g].nii_h_alpha     = d[g].nii_6548_h_alpha+d[g].nii_6584_h_alpha
          d[g].nii_h_alpha_err = sqrt(d[g].nii_6548_h_alpha_err^2+d[g].nii_6584_h_alpha_err^2)
       endif 

; ---------------------------------------------------------------------------
; [S II] 6716, 6731
       g = where((d.sii_6716_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].sii_6716_h_alpha = d[g].sii_6716_h_beta/ha[g]
          d[g].sii_6716_h_alpha_err = d[g].sii_6716_h_beta_err
;         d[g].sii_6716_h_alpha_err = im_compute_error(d[g].sii_6716_h_beta,$
;           d[g].sii_6716_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

       g = where((d.sii_6731_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].sii_6731_h_alpha = d[g].sii_6731_h_beta/ha[g]
          d[g].sii_6731_h_alpha_err = d[g].sii_6731_h_beta_err
;         d[g].sii_6731_h_alpha_err = im_compute_error(d[g].sii_6731_h_beta,$
;           d[g].sii_6731_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

       g = where((d.sii_6716_h_beta gt -900.0) and (d.sii_6731_h_beta gt -900.0),ng) ; compute the sum
       if (ng ne 0L) then begin
          d[g].sii_h_beta      = d[g].sii_6716_h_beta + d[g].sii_6731_h_beta
          d[g].sii_h_beta_err  = sqrt(d[g].sii_6716_h_beta_err^2.0 + d[g].sii_6731_h_beta_err^2.0)
       endif

       g = where((d.sii_6716_h_alpha gt -900.0) and (d.sii_6731_h_alpha gt -900.0),ng) ; compute the sum
       if (ng ne 0L) then begin
          d[g].sii_h_alpha      = d[g].sii_6716_h_alpha+d[g].sii_6731_h_alpha
          d[g].sii_h_alpha_err  = sqrt(d[g].sii_6716_h_alpha_err^2.0+d[g].sii_6731_h_alpha_err^2.0)
       endif

; ---------------------------------------------------------------------------
; [S III] 6312
       d.siii_6312_h_beta     = hiidata.siii_6312
       d.siii_6312_h_beta_err = hiidata.siii_6312_err
       g = where((d.siii_6312_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].siii_6312_h_alpha     = d[g].siii_6312_h_beta/ha[g]
          d[g].siii_6312_h_alpha_err = d[g].siii_6312_h_beta_err
;         d[g].siii_6312_h_alpha_err = im_compute_error(d[g].siii_6312_h_beta,$
;           d[g].siii_6312_h_beta_err,ha[g],ha_err[g],/quotient)
       endif

; ---------------------------------------------------------------------------
; [S III] 9069,9532

; from 9069 and 9532 separately
       g = where((hiidata.siii_9069 gt -900.0) and (hiidata.siii_9532 gt -900.0) and $
         (hiidata.siii_9069_9532 lt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].siii_9069_h_beta     = hiidata[g].siii_9069
          d[g].siii_9069_h_beta_err = hiidata[g].siii_9069_err
          d[g].siii_9532_h_beta     = hiidata[g].siii_9532
          d[g].siii_9532_h_beta_err = hiidata[g].siii_9532_err
          d[g].siii_h_beta          = hiidata[g].siii_9069+hiidata[g].siii_9532
          d[g].siii_h_beta_err      = sqrt(hiidata[g].siii_9069_err^2+hiidata[g].siii_9532_err^2)
       endif

; from just 9069       
       g = where((hiidata.siii_9069 gt -900.0) and (hiidata.siii_9532 lt -900.0) and $
         (hiidata.siii_9069_9532 lt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].siii_9069_h_beta     = hiidata[g].siii_9069
          d[g].siii_9069_h_beta_err = hiidata[g].siii_9069_err
          d[g].siii_9532_h_beta     = sratio*hiidata[g].siii_9069
          d[g].siii_9532_h_beta_err = sratio*hiidata[g].siii_9069_err
          d[g].siii_h_beta          = (1.0+sratio)*hiidata[g].siii_9069
          d[g].siii_h_beta_err      = (1.0+sratio)*hiidata[g].siii_9069_err
       endif

; from just 9532
       g = where((hiidata.siii_9069 lt -900.0) and (hiidata.siii_9532 gt -900.0) and $
         (hiidata.siii_9069_9532 lt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].siii_9069_h_beta     = (1.0/sratio)*hiidata[g].siii_9532
          d[g].siii_9069_h_beta_err = (1.0/sratio)*hiidata[g].siii_9532_err
          d[g].siii_9532_h_beta     = hiidata[g].siii_9532
          d[g].siii_9532_h_beta_err = hiidata[g].siii_9532_err
          d[g].siii_h_beta          = scor*hiidata[g].siii_9532
          d[g].siii_h_beta_err      = scor*hiidata[g].siii_9532_err
       endif

; from 9069 and 9532 together
       g = where((hiidata.siii_9069 lt -900.0) and (hiidata.siii_9532 lt -900.0) and $
         (hiidata.siii_9069_9532 gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].siii_h_beta          = hiidata[g].siii_9069_9532
          d[g].siii_h_beta_err      = hiidata[g].siii_9069_9532_err
          d[g].siii_9069_h_beta     = hiidata[g].siii_9069_9532/(1.0+sratio)
          d[g].siii_9069_h_beta_err = hiidata[g].siii_9069_9532_err/(1.0+sratio)/sqrt(2.0) ; note factor of 1.4
          d[g].siii_9532_h_beta     = hiidata[g].siii_9069_9532/scor
          d[g].siii_9532_h_beta_err = hiidata[g].siii_9069_9532_err/scor/sqrt(2.0) ; note factor of 1.4
       endif

       g = where((d.siii_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].siii_h_alpha     = d[g].siii_h_beta/ha[g]
          d[g].siii_h_alpha_err = d[g].siii_h_beta_err
;         d[g].siii_h_alpha_err = im_compute_error(d[g].siii_h_beta,$
;           d[g].siii_h_beta_err,ha[g],ha_err[g],/quotient)
       endif 

       g = where((d.siii_9069_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].siii_9069_h_alpha     = d[g].siii_9069_h_beta/ha[g]
          d[g].siii_9069_h_alpha_err = d[g].siii_9069_h_beta_err
;         d[g].siii_9069_h_alpha_err = im_compute_error(d[g].siii_9069_h_beta,$
;           d[g].siii_9069_h_beta_err,ha[g],ha_err[g],/quotient)
       endif 

       g = where((d.siii_9532_h_beta gt -900.0) and ha_good,ng)
       if (ng ne 0L) then begin
          d[g].siii_9532_h_alpha     = d[g].siii_9532_h_beta/ha[g]
          d[g].siii_9532_h_alpha_err = d[g].siii_9532_h_beta_err
;         d[g].siii_9532_h_alpha_err = im_compute_error(d[g].siii_9532_h_beta,$
;           d[g].siii_9532_h_beta_err,ha[g],ha_err[g],/quotient)
       endif 

; ---------------------------------------------------------------------------
; compute additional useful line ratios

; [N II] 6584 / [O II]       

       g = where((d.nii_6584_h_beta gt -900.0) and (d.oii_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].nii_6584_oii     = d[g].nii_6584_h_beta/d[g].oii_h_beta
          d[g].nii_6584_oii_err = im_compute_error(d[g].nii_6584_h_beta,$
            d[g].nii_6584_h_beta_err,d[g].oii_h_beta,d[g].oii_h_beta_err,/quotient)
       endif

; [N II] 6584 / [O III] 5007

       g = where((d.nii_6584_h_beta gt -900.0) and (d.oiii_5007_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].nii_6584_oiii_5007     = d[g].nii_6584_h_beta/d[g].oiii_5007_h_beta
          d[g].nii_6584_oiii_5007_err = im_compute_error(d[g].nii_6584_h_beta,$
            d[g].nii_6584_h_beta_err,d[g].oiii_5007_h_beta,d[g].oiii_5007_h_beta_err,/quotient)
       endif
       
; [N II] 6584 / [S II] 6716+6731

       g = where((d.nii_6584_h_beta gt -900.0) and (d.sii_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].nii_6584_sii     = d[g].nii_6584_h_beta/d[g].sii_h_beta
          d[g].nii_6584_sii_err = im_compute_error(d[g].nii_6584_h_beta,$
            d[g].nii_6584_h_beta_err,d[g].sii_h_beta,d[g].sii_h_beta_err,/quotient)
       endif
              
; [O II] / [N II] 6584

       g = where((d.oii_h_beta gt -900.0) and (d.nii_6584_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oii_nii_6584     = d[g].oii_h_beta/d[g].nii_6584_h_beta
          d[g].oii_nii_6584_err = im_compute_error(d[g].oii_h_beta,d[g].oii_h_beta_err,$
            d[g].nii_6584_h_beta,d[g].nii_6584_h_beta_err,/quotient)
       endif

; [O II] / [O III] 5007

       g = where((d.oii_h_beta gt -900.0) and (d.oiii_5007_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oii_oiii_5007     = d[g].oii_h_beta/d[g].oiii_5007_h_beta
          d[g].oii_oiii_5007_err = im_compute_error(d[g].oii_h_beta,d[g].oii_h_beta_err,$
            d[g].oiii_5007_h_beta,d[g].oiii_5007_h_beta_err,/quotient)
       endif

; [O II] / [S II] 6716+6731
       
       g = where((d.oii_h_beta gt -900.0) and (d.sii_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oii_sii     = d[g].oii_h_beta/d[g].sii_h_beta
          d[g].oii_sii_err = im_compute_error(d[g].oii_h_beta,d[g].oii_h_beta_err,$
            d[g].sii_h_beta,d[g].sii_h_beta_err,/quotient)
       endif

; ([O II]/Hb) * ([O II]/[S II] 6716+6731)
       
       g = where((d.oii_h_beta gt -900.0) and (d.oii_sii gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oii_h_beta_oii_sii = d[g].oii_h_beta*d[g].oii_sii
          d[g].oii_h_beta_oii_sii_err = im_compute_error(d[g].oii_h_beta,$
            d[g].oii_h_beta_err,d[g].oii_sii,d[g].oii_sii_err,/product)
       endif

; [O III] 5007 / [N II] 6584
       
       g = where((d.oiii_5007_h_beta gt -900.0) and (d.nii_6584_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oiii_5007_nii_6584     = d[g].oiii_5007_h_beta/d[g].nii_6584_h_beta
          d[g].oiii_5007_nii_6584_err = im_compute_error(d[g].oiii_5007_h_beta,$
            d[g].oiii_5007_h_beta_err,d[g].nii_6584_h_beta,d[g].nii_6584_h_beta_err,/quotient)
       endif

; [O III] 5007 / [O II]
       
       g = where((d.oiii_5007_h_beta gt -900.0) and (d.oii_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oiii_5007_oii     = d[g].oiii_5007_h_beta/d[g].oii_h_beta
          d[g].oiii_5007_oii_err = im_compute_error(d[g].oiii_5007_h_beta,$
            d[g].oiii_5007_h_beta_err,d[g].oii_h_beta,d[g].oii_h_beta_err,/quotient)
       endif
       
; [O III] 5007 / [S II] 6716+6731
       
       g = where((d.oiii_5007_h_beta gt -900.0) and (d.sii_h_beta gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oiii_5007_sii     = d[g].oiii_5007_h_beta/d[g].sii_h_beta
          d[g].oiii_5007_sii_err = im_compute_error(d[g].oiii_5007_h_beta,$
            d[g].oiii_5007_h_beta_err,d[g].sii_h_beta,d[g].sii_h_beta_err,/quotient)
       endif

; ([O III] 5007 / Hb) / ([N II] 6584 / Ha)
       
       g = where((d.oiii_5007_h_beta gt -900.0) and (d.nii_6584_h_alpha gt -900.0),ng)
       if (ng ne 0L) then begin
          d[g].oiii_5007_h_beta_nii_6584_h_alpha     = d[g].oiii_5007_h_beta/d[g].nii_6584_h_alpha
          d[g].oiii_5007_h_beta_nii_6584_h_alpha_err = im_compute_error(d[g].oiii_5007_h_beta,$
            d[g].oiii_5007_h_beta_err,d[g].nii_6584_h_alpha,d[g].nii_6584_h_alpha_err,/quotient)
       endif

; O32, R23, and S23 get computed in WRITE_HII_REGIONS by
; IM_ABUNDANCE(), so don't bother computing them here       
       
;; O32 = [O III] 4959, 5007 / [O II] 3727
;
;       g = where((d.oii_h_beta gt -900.0) and (d.oiii_h_beta gt -900.0),ng)
;       if (ng ne 0L) then begin
;          d[g].o32 = d[g].oiii_h_beta/d[g].oii_h_beta
;          d[g].o32_err = im_compute_error(d[g].oiii_h_beta,d[g].oiii_h_beta_err,$
;            d[g].oii_h_beta,d[g].oii_h_beta_err,/quotient)
;       endif
;
;; R23 = [O II] 3727 + [O III] 4959, 5007
;
;       g = where((d.oii_h_beta gt -900.0) and (d.oiii_h_beta gt -900.0),ng)
;       if (ng ne 0L) then begin
;          d[g].r23 = d[g].oii_h_beta + d[g].oiii_h_beta
;          d[g].r23_err = sqrt(d[g].oii_h_beta_err^2.0 + d[g].oiii_h_beta_err^2.0)
;       endif
;
;; S23 = [S II] 6716,6731 + [S III] 9069,9532
;
;       g = where((d.sii_h_beta gt -900.0) and (d.siii_h_beta gt -900.0),ng)
;       if (ng ne 0L) then begin
;          d[g].s23 = d[g].sii_h_beta + d[g].siii_h_beta
;          d[g].s23_err = sqrt(d[g].sii_h_beta_err^2.0 + d[g].siii_h_beta_err^2.0)
;       endif

; ---------------------------------------------------------------------------       
; before taking the logarithm of each line flux, initialize the
; linefit data structure
       linefit1 = init_hiiregions_linefit(d)

; ---------------------------------------------------------------------------       
; compute the logarithm of all the line fluxes; also return the data
; structure where no logarithms have been taken
       dnolog = d

       except = ['HII_*','NED_*','GALAXY_*','EWHB',$
         'FHA','TEXREF','REFERENCE','*_ERR','LIT_*']
       sub = struct_trimtags(d,except=except)

       tags = tag_names(d)
       linetags = tag_names(sub)
       
       ntags = n_elements(linetags)
       for k = 0L, ntags-1L do begin

          match1 = where(strmatch(tags,linetags[k]) eq 1B)
          match2 = where(strmatch(tags,linetags[k]+'_ERR') eq 1B)
          
          flux = d.(match1)
          ferr = d.(match2)

          newflux = flux
          newferr = ferr
          
          good = where(flux gt -900.0,ngood)
          if (ngood ne 0L) then begin
;            print, k, minmax(flux[good])
             zero = where(flux[good] eq 0.0,nzero)
             if (nzero ne 0L) then message, 'Problem!'
             newferr[good] = ferr[good]/flux[good]/alog(10.0)
             newflux[good] = alog10(flux[good])
          endif

          d.(match1) = newflux
          d.(match2) = newferr
       endfor

       if (j eq 0L) then begin
          data = d 
          linefit = linefit1
          data_nolog = dnolog
       endif else begin
          data = struct_append(data,d)
          linefit = struct_append(linefit,linefit1)
          data_nolog = struct_append(data_nolog,dnolog)
       endelse

       print, file_basename(datafile[j])
;      print, file_basename(datafile[j]) & stop ; cc = get_kbrd(1)

    endfor 

; done!    
    data = reform(data)
    linefit = reform(linefit)
    data_nolog = reform(data_nolog)
       
return, data
end
