; ###########################################################################
; THIS HAS LOTS OF CALCULATIONS THAT MIGHT BE USEFUL
; ###########################################################################
;+
; NAME:
;       IUNRED_LINEDUST()
;
; PURPOSE:
;       Correct emission-line fluxes for dust extinction using the
;       observed Balmer decrement.
;
; INPUTS:
;       linedust - input data structure
;
; OPTIONAL INPUTS:
;       ratio           - E(B-V)_continuum divided by E(B-V)_gas
;                         (default 1.0, in general should be less than
;                         1.0) 
;       snrcut_linedust - 
;       cutsig          - 
;       extra           - keywords for GET_EBV()
;
; KEYWORD PARAMETERS:
;       HaHb        - de-redden using the H-alpha/H-beta line ratio 
;                     (default) 
;       HbHg        - de-redden using the H-beta/H-gamma line ratio 
;       HaHg        - de-redden using the H-alpha/H-gamma line ratio 
;       combination - correct for extinction using either the Ha/Hb or
;                     the Hb/Hg ratio, giving preference to Ha/Hb
;       nopropagate - do not propagate the error in the reddening
;                     correction to the unreddened fluxes
;       allow_r_min - use the minimum physically plausible Balmer
;                     decrement, rather than the mean (T=10^4 K) value 
;       silent      - suppress messages to STDOUT
;
; OUTPUTS:
;       linedust   - input data structure (modified)
;       linenodust - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Upper limits on the Balmer lines are not treated properly. 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 November 6, U of A
;       jm03jul23uofa - added EHAHB and EHAHG structure fields 
;       jm03aug3uofa  - properly correct the emission line luminosities
;       jm03aug6uofa  - added temperature dependence on the Balmer
;                       decrement and HBETA keyword
;       jm04jan13uofa - added SNRCUT_LINEDUST keyword and additional error
;                       checking 
;       jm04jul29uofa - additional error checking
;       jm04sep20uofa - added NOPROPAGATE keyword
;       jm04nov01uofa - correct the absorption-uncorrected Balmer
;                       emission-line fluxes for reddening (see
;                       PARSE_ISPECLINEFIT) 
;       jm04nov02uofa - added the H_GAMMA_PREDICT and H_DELTA_PREDICT
;                       structures to LINEDUST
;       jm04nov03uofa - removed HGAMMA and HBETA keywords in favor of
;                       HaHb, HaHg, and HbHg keywords; use new tag
;                       names for the reddening-corrected EW's 
;       jm05feb02uofa - removed TEMPERATURE keyword; this keyword is
;                       now controlled entirely in GET_EBV() through
;                       the _EXTRA keyword
;       jm05nov13uofa - hard-wire the Balmer decrements
;       jm06apr12uofa - CUTSIG optional input added
;       jm08feb06nyu  - added ALLOW_R_MIN keyword
;       jm08feb20nyu  - this routine was in desperate need of being
;                       streamlined, so for now remove the Ha/Hg
;                       calculations
;
; Copyright (C) 2002-2006, 2008, John Moustakas
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

function iunred_linedust, linedust, ratio=ratio, snrcut_linedust=snrcut_linedust, $
  cutsig=cutsig, HaHb=HaHb, HbHg=HbHg, HaHg=HaHg, combination=combination, $
  nopropagate=nopropagate, allow_r_min=allow_r_min, silent=silent, _extra=extra
    
    lsun = 3.826D33  ; [erg/s]

    nspec = n_elements(linedust)
    if (nspec eq 0L) then begin
       doc_library, 'iunred_linedust'
       return, -1L
    endif

    if (n_elements(ratio) eq 0L) then ratio = 1.0
    if (n_elements(snrcut_linedust) eq 0L) then snrcut_linedust = 5.0
    if (n_elements(cutsig) eq 0L) then cutsig = 1.0

    if not keyword_set(silent) then splog, 'S/N > '+string(snrcut_linedust,format='(G0.0)')+'.'

    if keyword_set(hahg) then begin
       splog, 'HAHG keyword temporarily relegated!'
       return, linedust
    endif
    
    temp_max = 20000.0 ; maximum temperature
    
    linename = strlowcase(strcompress(linedust[0].linename,/remove))
    nline = n_elements(linename)

    R_HaHb = 2.86 ; intrinsic Ha/Hb ratio
    R_HbHg = 2.14 ; intrinsic Hb/Hg ratio
    R_HaHg = 6.11 ; intrinsic Ha/Hg ratio
    R_HaHd = 11.1 ; intrinsic Ha/Hd ratio

;   R_HaHb = return_tbalmer(_extra,/HaHb) ; intrinsic Ha/Hb ratio
;   R_HbHg = return_tbalmer(_extra,/HbHg) ; intrinsic Hb/Hg ratio
;   R_HaHg = return_tbalmer(_extra,/HaHg) ; intrinsic Ha/Hg ratio
;   R_HaHd = return_tbalmer(_extra,/HaHd) ; intrinsic Ha/Hd ratio

; minimum physically plausible Balmer decrements
    
    R_HaHb_min = 2.75
    R_HbHg_min = 2.10
    R_HaHg_min = 5.78
    R_HaHd_min = 10.4
    
;   R_HaHb_min = return_tbalmer(temp_max,/HaHb)
;   R_HbHg_min = return_tbalmer(temp_max,/HbHg)
;   R_HaHg_min = return_tbalmer(temp_max,/HaHg)
;   R_HaHd_min = return_tbalmer(temp_max,/HaHd)
    
; initialize and fill the output structure

    newtags = {$
      hahb:  0.0, hahb_err:  -1.0, ebv_hahb: 0.0, ebv_hahb_err: -1.0, ebv_hahb_negative: 0B, $
      hbhg:  0.0, hbhg_err:  -1.0, ebv_hbhg: 0.0, ebv_hbhg_err: -1.0, ebv_hbhg_negative: 0B, $
;     hahg:  0.0, hahg_err:  -1.0, ebv_hahg: 0.0, ebv_hahg_err: -1.0, ebv_hahg_negative: 0B, $
      ehbha: 0.0, ehbha_err: -1.0, ehghb:    0.0, ehghb_err:    -1.0}
;     ehgha: 0.0, ehgha_err: -1.0}
    newtags = replicate(newtags,nspec)
      
    linenodust = struct_addtags(linedust,newtags)

; ---------------------------------------------------------------------------
; Ha/Hb decrement, corresponding E(B-V), and E(Hb-Ha) 
; ---------------------------------------------------------------------------

    if tag_exist(linedust,'H_ALPHA') and tag_exist(linedust,'H_BETA') then begin
       
; keep all objects with well-measured H-alpha and H-beta lines

       goodsnr = where((linenodust.h_alpha[1] gt 0.0) and $
         (linenodust.h_alpha[0]/linenodust.h_alpha[1] gt snrcut_linedust) and $
         (linenodust.h_beta[1] gt 0.0) and $
         (linenodust.h_beta[0]/linenodust.h_beta[1] gt snrcut_linedust),ngood)

       if (ngood ne 0L) then begin

; compute the Balmer decrement          
          
          linenodust[goodsnr].hahb = linenodust[goodsnr].h_alpha[0]/linenodust[goodsnr].h_beta[0]
          linenodust[goodsnr].hahb_err = im_compute_error(linenodust[goodsnr].h_alpha[0],$
            linenodust[goodsnr].h_alpha[1],linenodust[goodsnr].h_beta[0],linenodust[goodsnr].h_beta[1],$
            /quotient)

;         niceprint, linenodust[goodsnr].hahb, linenodust[goodsnr].hahb_err

; do not compute E(B-V) for objects with unphysically low Balmer
; decrements; they must be within CUTSIG sigma of the adopted/minimum
; decrement; there is no upper limit on the Balmer decrement (dusty!);
; note, however, that objects can have negative E(Hb-Ha) that is zero
; within the errors; however, for these objects keep the measured
; Balmer decrements (ie, in general they will be less than 2.86)

; jm04jul22uofa - for objects with E(B-V)<0 set E(B-V)=0 and flag them
          
; jm05feb03uofa - new conservative policy: require the decrement to be
;                 within CUTSIG of the **adopted** intrinsic
;                 decrement, **not** the minimum possible decrement

; jm08feb06nyu - use a keyword to allow the smallest possible decrement
          
          if keyword_set(allow_r_min) then $
            keep = where(linenodust[goodsnr].hahb+cutsig*linenodust[goodsnr].hahb_err ge R_HaHb_min,nkeep) else $
              keep = where(linenodust[goodsnr].hahb+cutsig*linenodust[goodsnr].hahb_err ge R_HaHb,nkeep)

          if (nkeep ne 0L) then begin

;            srt = reverse(sort(linenodust[goodsnr[keep]].hahb))
;            niceprint, linenodust[goodsnr[keep][srt]].galaxy, linenodust[goodsnr[keep][srt]].hahb, $
;              linenodust[goodsnr[keep][srt]].hahb_err, linenodust[goodsnr[keep][srt]].hahb/R_HaHb
             
             linenodust[goodsnr[keep]].ebv_hahb = get_ebv(linenodust[goodsnr[keep]].hahb,$
               decrement_err=linenodust[goodsnr[keep]].hahb_err,ebv_err=ebv_err,color=color,$
               err_color=color_err,/HaHb,_extra=extra)

             linenodust[goodsnr[keep]].ebv_hahb_err = ebv_err
             linenodust[goodsnr[keep]].ehbha = color
             linenodust[goodsnr[keep]].ehbha_err = color_err

          endif

          neg = where(linenodust.ebv_hahb lt 0.0,nneg)
          good = where(linenodust.ebv_hahb_err gt 0.0,ngood)

          if (nneg ne 0L) then begin
             if (not keyword_set(silent)) then splog, 'Identified '+string(nneg,format='(I0)')+'/'+$
               string(ngood,format='(I0)')+' [of '+string(nspec,format='(I0)')+$
               ' objects] with negative E(B-V) [Ha/Hb].'
             linenodust[neg].ebv_hahb = 0.0
             linenodust[neg].ebv_hahb_negative = 1B
          endif 
          
       endif else splog, 'WARNING: No good Ha/Hb Balmer decrement measurements!'

    endif

; ---------------------------------------------------------------------------
; Hb/Hg decrement, corresponding E(B-V), and E(Hg-Hb) 
; ---------------------------------------------------------------------------

    if tag_exist(linedust,'H_GAMMA') and tag_exist(linedust,'H_BETA') then begin
    
; keep all objects with well-measured H-gamma and H-beta lines
       
       goodsnr = where((linenodust.h_gamma[1] gt 0.0) and $
         (linenodust.h_gamma[0]/linenodust.h_gamma[1] gt snrcut_linedust) and $
         (linenodust.h_beta[1] gt 0.0) and $
         (linenodust.h_beta[0]/linenodust.h_beta[1] gt snrcut_linedust),ngood)

       if (ngood ne 0L) then begin

; compute the Balmer decrement          
          
          linenodust[goodsnr].hbhg = linenodust[goodsnr].h_beta[0]/linenodust[goodsnr].h_gamma[0]
          linenodust[goodsnr].hbhg_err = im_compute_error(linenodust[goodsnr].h_beta[0],$
            linenodust[goodsnr].h_beta[1],linenodust[goodsnr].h_gamma[0],linenodust[goodsnr].h_gamma[1],$
            /quotient)

; only keep objects with physical decrements; see notes from 05feb03 above

          if keyword_set(allow_r_min) then $
            keep = where(linenodust[goodsnr].hbhg+cutsig*linenodust[goodsnr].hbhg_err ge R_HbHg_min,nkeep) else $
              keep = where(linenodust[goodsnr].hbhg+cutsig*linenodust[goodsnr].hbhg_err ge R_HbHg,nkeep)

          if (nkeep ne 0L) then begin

             linenodust[goodsnr[keep]].ebv_hbhg = get_ebv(linenodust[goodsnr[keep]].hbhg,$
               decrement_err=linenodust[goodsnr[keep]].hbhg_err,ebv_err=ebv_err,color=color,$
               err_color=color_err,/HbHg,_extra=extra)

             linenodust[goodsnr[keep]].ebv_hbhg_err = ebv_err
             linenodust[goodsnr[keep]].ehghb = color
             linenodust[goodsnr[keep]].ehghb_err = color_err

          endif
             
          neg = where(linenodust.ebv_hbhg lt 0.0,nneg)
          good = where(linenodust.ebv_hbhg_err gt 0.0,ngood)

          if (nneg ne 0L) then begin
             if (not keyword_set(silent)) then splog, 'Identified '+string(nneg,format='(I0)')+'/'+$
               string(ngood,format='(I0)')+' [of '+string(nspec,format='(I0)')+$
               ' objects] with negative E(B-V) [Hb/Hg].'
             linenodust[neg].ebv_hbhg = 0.0
             linenodust[neg].ebv_hbhg_negative = 1B
          endif 
          
       endif else splog, 'WARNING: No good Hb/Hg Balmer decrement measurements!'
          
    endif

; #########################    
; THIS CODE IS GOOD!!
; #########################    
;; ---------------------------------------------------------------------------
;; Ha/Hg decrement, corresponding E(B-V), and E(Hg-Ha) 
;; ---------------------------------------------------------------------------
;
;    if tag_exist(linedust,'H_GAMMA') and tag_exist(linedust,'H_ALPHA') then begin
;    
;; keep all objects with well-measured H-gamma and H-alpha lines
;       
;       goodsnr = where((linenodust.h_gamma[1] gt 0.0) and $
;         (linenodust.h_gamma[0]/linenodust.h_gamma[1] gt snrcut_linedust) and $
;         (linenodust.h_alpha[1] gt 0.0) and $
;         (linenodust.h_alpha[0]/linenodust.h_alpha[1] gt snrcut_linedust),ngood)
;
;       if (ngood ne 0L) then begin
;
;; compute the Balmer decrement          
;          
;          linenodust[goodsnr].hahg = linenodust[goodsnr].h_alpha[0]/linenodust[goodsnr].h_gamma[0]
;          linenodust[goodsnr].hahg_err = im_compute_error(linenodust[goodsnr].h_alpha[0],$
;            linenodust[goodsnr].h_alpha[1],linenodust[goodsnr].h_gamma[0],linenodust[goodsnr].h_gamma[1],$
;            /quotient)
;
;; only keep objects with physical decrements; see notes from 05feb03 above
;
;          if keyword_set(allow_r_min) then $
;            keep = where(linenodust[goodsnr].hahg+cutsig*linenodust[goodsnr].hahg_err ge R_HaHg_min,nkeep) else $
;              keep = where(linenodust[goodsnr].hahg+cutsig*linenodust[goodsnr].hahg_err ge R_HaHg,nkeep)
;
;          if (nkeep ne 0L) then begin
;
;             linenodust[goodsnr[keep]].ebv_hahg = get_ebv(linenodust[goodsnr[keep]].hahg,$
;               decrement_err=linenodust[goodsnr[keep]].hahg_err,ebv_err=ebv_err,color=color,$
;               err_color=color_err,/HaHg,_extra=extra)
;
;             linenodust[goodsnr[keep]].ebv_hahg_err = ebv_err
;             linenodust[goodsnr[keep]].ehgha = color
;             linenodust[goodsnr[keep]].ehgha_err = color_err
;
;          endif
;             
;          neg = where(linenodust.ebv_hahg lt 0.0,nneg)
;          good = where(linenodust.ebv_hahg_err gt 0.0,ngood)
;
;          if (nneg ne 0L) then begin
;             if (not keyword_set(silent)) then splog, 'Identified '+string(nneg,format='(I0)')+'/'+$
;               string(ngood,format='(I0)')+' [of '+string(nspec,format='(I0)')+$
;               ' objects] with negative E(B-V) [Ha/Hg].'
;             linenodust[neg].ebv_hahg = 0.0
;             linenodust[neg].ebv_hahg_negative = 1B
;          endif
;          
;       endif else splog, 'WARNING: No good Ha/Hg Balmer decrement measurements!'
;          
;    endif
       
; ---------------------------------------------------------------------------

; correct the emission-line fluxes for dust extinction; if COMBINATION
; is set then use Hb/Hg where Ha/Hb is not defined (e.g., for
; higher-redshift objects)

    if keyword_set(combination) then begin

       good = where((linenodust.ebv_hahb_err gt -1.0) or $
         (linenodust.ebv_hbhg_err gt -1.0),ngood,comp=bad,ncomp=nbad)

       if (ngood ne 0L) then begin

          ebv = dblarr(ngood)
          ebv_err = dblarr(ngood)

          ha = where((linenodust[good].ebv_hahb_err gt -1.0),nha)
          hb = where((linenodust[good].ebv_hbhg_err gt -1.0),nhb)

          if (nha ne 0L) then begin
             ebv[ha] = linenodust[good[ha]].ebv_hahb
             ebv_err[ha] = linenodust[good[ha]].ebv_hahb_err
          endif

          if (nhb ne 0L) then begin
             ebv[hb] = linenodust[good[hb]].ebv_hbhg
             ebv_err[hb] = linenodust[good[hb]].ebv_hbhg_err
          endif

       endif else begin

          splog, 'Ha/Hb, Hb/Hg: No good Balmer decrements!'
          return, linenodust

       endelse

    endif

    if keyword_set(HbHg) then begin

       good = where(linenodust.ebv_hbhg_err gt -1.0,ngood,comp=bad,ncomp=nbad)

       if (ngood ne 0L) then begin
          ebv = linenodust[good].ebv_hbhg
          ebv_err = linenodust[good].ebv_hbhg_err
       endif else begin
          splog, 'Hb/Hg: No good Balmer decrements!'
          return, linenodust
       endelse

    endif

; #########################    
; THIS CODE IS GOOD!!
; #########################    
;   if keyword_set(HaHg) then begin
;
;      good = where(linenodust.ebv_hahg_err gt -1.0,ngood,comp=bad,ncomp=nbad)
;
;      if (ngood ne 0L) then begin
;         ebv = linenodust[good].ebv_hahg
;         ebv_err = linenodust[good].ebv_hahg_err
;      endif else begin
;         splog, 'Ha/Hg: No good Balmer decrements!'
;         return, linenodust
;      endelse
;
;   endif

    if (not keyword_set(combination)) and (not keyword_set(HbHg)) and $
      (not keyword_set(HaHg)) then HaHb = 1L

    if keyword_set(HaHb) then begin 

       good = where(linenodust.ebv_hahb_err gt -1.0,ngood,comp=bad,ncomp=nbad)

       if (ngood ne 0L) then begin
          ebv = linenodust[good].ebv_hahb
          ebv_err = linenodust[good].ebv_hahb_err
       endif else begin
          splog, 'Ha/Hb: No good Balmer decrements!'
          return, linenodust
       endelse

    endif

    if keyword_set(nopropagate) then ebv_err = ebv*0.0 ; no error!

; use different tag names for the reddening-corrected equivalent
; widths; we do not touch the _UNCOR and _EW_UNCOR tags created by
; PARSE_ISPECLINEFIT() 

    newtags = create_struct(linename[0]+'_ew_cor',[0.0,-2.0])
    for iline = 1L, nline-1L do newtags = create_struct(newtags,linename[iline]+'_ew_cor', [0.0,-2.0])
    newtags = replicate(newtags,nspec)

    linenodust = struct_addtags(temporary(linenodust),newtags)
    
; the T structure data points can be corrected for reddening but the
; NORED structure data point errors must be set to -1.0
    
    t = linenodust[good]
    if (nbad ne 0L) then nored = linenodust[bad]

; emission-line luminosity fields

    if tag_exist(t,'LINELUMNAME') then begin

       lumtags = 1 ; true
       linelumname = strlowcase(strcompress(t[0].linelumname,/remove))
       nlinelum = n_elements(linelumname)

    endif else lumtags = 0 ; false

; loop on each line and correct for dust.  also, if emission line
; luminosities have been computed then correct those for dust as well;
; also correct the emission-line continuum and form new
; "dust-corrected" equivalent widths (in the general case where RATIO
; is not equal to 1.0)

    tags = tag_names(linenodust)
    for i = 0L, nline-1L do begin

; only de-redden well-measured emission lines.  should the upper
; limits be de-reddened as well?  don't for now.

       match = where(linename[i] eq strlowcase(tags),nmatch)
       if (nmatch eq 0L) then begin
          splog, 'Required structure tag '+strupcase(linename[i])+' not found.'
          return, -1L
       endif

       lineflux = reform(double(t.(match)[0,*]))
       lineflux_err = reform(double(t.(match)[1,*]))

       wavematch = where(linename[i]+'_wave' eq strlowcase(tags),nwavematch)
       if (nwavematch eq 0L) then begin
          splog, 'Required structure tag '+strupcase(linename[i])+'_WAVE not found.'
          return, -1L
       endif

       linewave = t.(wavematch)

       nice = where(lineflux_err gt 0.0,nnice,comp=poor,ncomp=npoor)

       newflux = lineflux
       newflux_err = lineflux_err

       if (nnice gt 0L) then begin

          newflux[nice] = dust_correct(lineflux[nice],linewave[nice],$
            err_lineflux=lineflux_err[nice],ebv=ebv[nice],err_ebv=ebv_err[nice],$
            err_newflux=tempflux_err)
          newflux_err[nice] = tempflux_err

       endif

       nodust = t.(match)
       nodust = transpose([ [newflux], [newflux_err] ])
;      nodust[*,nice] = transpose([ [newflux], [newflux_err] ])
       t.(match) = float(nodust)
       if (nbad ne 0L) then nored.(match) = [0.0,-1.0]

; next de-redden the continuum by an amount E(B-V)_continuum that can
; differ from E(B-V)_gas by a factor RATIO

       match = where(linename[i]+'_continuum' eq strlowcase(tags),nmatch)
       if (nmatch ne 0L) then begin
          
          contflux = reform(double(t.(match)[0,*]))
          contflux_err = reform(double(t.(match)[1,*]))

          nice = where(contflux_err gt 0.0,nnice,comp=poor,ncomp=npoor)

          newcontflux = contflux
          newcontflux_err = contflux_err

          ebv_cont = ebv*ratio

          if (nnice gt 0L) then begin

             newcontflux[nice] = dust_correct(contflux[nice],linewave[nice],$
               err_lineflux=contflux_err[nice],ebv=ebv_cont[nice],err_ebv=ebv_err[nice],$
               err_newflux=tempflux_err)
             newcontflux_err[nice] = tempflux_err

          endif

          nodust = t.(match)
          nodust = transpose([ [newcontflux], [newcontflux_err] ])
          t.(match) = float(nodust)
          if (nbad ne 0L) then nored.(match) = [0.0,-1.0]

       endif else splog, 'Structure tag '+strupcase(linename[i])+'_CONTINUUM not present.'
       
; finally re-compute the equivalent widths

       match = where(linename[i]+'_ew' eq strlowcase(tags),nmatch)
       cormatch = where(linename[i]+'_ew_cor' eq strlowcase(tags),ncormatch)
       if (nmatch ne 0L) then begin

          ew = reform(double(t.(match)[0,*]))
          ew_err = reform(double(t.(match)[1,*]))
          
          nice = where((newflux_err gt 0.0) and (newcontflux_err gt 0.0),nnice)

          new_ew = ew
          new_ew_err = ew_err

          if (nnice ne 0L) then begin

             new_ew[nice] = newflux[nice]/newcontflux[nice]
             new_ew_err[nice] = im_compute_error(newflux[nice],newflux_err[nice],$
               newcontflux[nice],newcontflux_err[nice],/quotient)

          endif
          
          nodust = t.(cormatch)
          nodust = transpose([ [new_ew], [new_ew_err] ])
          t.(cormatch) = float(nodust)
          if (nbad ne 0L) then nored.(cormatch) = [0.0,-1.0]

       endif else splog, 'Structure tag '+strupcase(linename[i])+'_EW not present.'

; if this line is a Balmer line, and there exists an absorption
; uncorrected fluxes from PARSE_ISPECLINEFIT(), then correct the flux
; and the equivalent width for reddening

       uncorlinename = linename[i]+'_uncor'
       uncorlineewname = linename[i]+'_ew_uncor'

       if tag_exist(linedust,uncorlinename) then begin

          match = where(uncorlinename eq strlowcase(tags),nmatch)
          
          lineflux = reform(double(t.(match)[0,*]))
          lineflux_err = reform(double(t.(match)[1,*]))
          
          nice = where(lineflux_err gt 0.0,nnice,comp=poor,ncomp=npoor)

          newflux = lineflux
          newflux_err = lineflux_err

          if (nnice gt 0L) then begin

             newflux[nice] = dust_correct(lineflux[nice],linewave[nice],$
               err_lineflux=lineflux_err[nice],ebv=ebv[nice],err_ebv=ebv_err[nice],$
               err_newflux=tempflux_err)
             newflux_err[nice] = tempflux_err

          endif

          nodust = t.(match)
          nodust = transpose([ [newflux], [newflux_err] ])

          t.(match) = float(nodust)
          if (nbad ne 0L) then nored.(match) = [0.0,-1.0]

; re-compute the equivalent width

          match = where(uncorlineewname eq strlowcase(tags),nmatch)

          ew = reform(double(t.(match)[0,*]))
          ew_err = reform(double(t.(match)[1,*]))
          
          nice = where((newflux_err gt 0.0) and (newcontflux_err gt 0.0),nnice)

          new_ew = ew
          new_ew_err = ew_err

          if (nnice ne 0L) then begin

             new_ew[nice] = newflux[nice]/newcontflux[nice]
             new_ew_err[nice] = im_compute_error(newflux[nice],newflux_err[nice],$
               newcontflux[nice],newcontflux_err[nice],/quotient)

          endif
          
          nodust = t.(match)
          nodust = transpose([ [new_ew], [new_ew_err] ])
          t.(match) = float(nodust)
          if (nbad ne 0L) then nored.(match) = [0.0,-1.0]
          
       endif
       
   endfor 

; emission-line luminosities

    if lumtags then begin

       lums = struct_trimtags(t,select_tags=linelumname+'_LUM')

       for j = 0L, nlinelum-1L do begin
          
          match = where(linelumname[j]+'_lum' eq strlowcase(tags))
          linelum = reform(double(lums.(j)[0,*]))
          linelum_err = reform(double(lums.(j)[1,*]))

          wavematch = where(linelumname[j]+'_wave' eq strlowcase(tags))
          linewave = t.(wavematch)

; only de-redden well-measured emission lines.  should the upper
; limits be de-reddened as well?  don't for now.

          nice = where(linelum_err gt 0.0,nnice,comp=poor,ncomp=npoor)

          newlum = linelum
          newlum_err = linelum_err

          if (nnice gt 0L) then begin

             linelum[nice] = lsun*10.0^(linelum[nice])
             linelum_err[nice] = alog(10.0)*linelum[nice]*linelum_err[nice]
             
             newlum[nice] = dust_correct(linelum[nice],linewave[nice],$
               err_lineflux=linelum_err[nice],ebv=ebv[nice],err_ebv=ebv_err[nice],$
               err_newflux=templum_err)
             newlum_err[nice] = templum_err
             
             newlum_err[nice] = newlum_err[nice]/newlum[nice]/alog(10.0)
             newlum[nice] = alog10(newlum[nice]/lsun)

          endif

          nodust = lums.(j)
          nodust = transpose([ [newlum], [newlum_err] ])
;         nodust[*,nice] = transpose([ [newlum], [newlum_err] ])

          t.(match) = float(nodust)
          if (nbad ne 0L) then nored.(match) = [0.0,-1.0]

       endfor 
       
    endif 

; now add the line flux corrections to LINENODUST

    linenodust[good] = struct_trimtags(t,select_tags=tag_names(t))
    if (nbad ne 0L) then linenodust[bad] = struct_trimtags(nored,select_tags=tag_names(nored))

; predict the H-beta, H-gamma, and H-delta reddened and unreddened
; fluxes by scaling the H-alpha flux and assuming case B

    if tag_exist(linedust,'H_ALPHA') then begin

       good = where((linenodust.h_alpha[1] gt 0.0) and (linenodust.ebv_hahb_err gt 0.0),ngood)

; H-BETA

       kl = k_lambda(4860.80,_extra=extra)
       
       hb = {$
         H_BETA_DUST_PREDICT:   [0.0,-2.0], $
         H_BETA_NODUST_PREDICT: [0.0,-2.0], $
         EBV_HAHB_PREDICT:      0.0, $
         EBV_HAHB_ERR_PREDICT: -1.0}
       hb = replicate(hb,nspec)

       if (ngood ne 0L) then begin

          ebv = linenodust[good].ebv_hahb
          ebv_err = linenodust[good].ebv_hahb_err

          flux = linenodust[good].h_alpha[0] / R_HaHb
          ferr = linenodust[good].h_alpha[1] / R_HaHb
          
          hb[good].h_beta_nodust_predict[0] = flux
          hb[good].h_beta_nodust_predict[1] = ferr

          hb[good].h_beta_dust_predict[0] = flux * 10^(-0.4*ebv*kl)
          hb[good].h_beta_dust_predict[1] = sqrt( (ferr*10.0^(0.4*ebv*kl))^2.0 + $
            (flux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*ebv_err)^2.0 )

;         plot, alog10(hb[good].h_beta_dust_predict[0]), linedust[good].h_beta[0] / $
;           hb[good].h_beta_dust_predict[0], xsty=3, ysty=3, ps=4
;         djs_oplot, !x.crange, [1,1], line=0, thick=2.0

          hahb = linedust[good].h_alpha[0]/hb[good].h_beta_dust_predict[0]
          hahb_err = im_compute_error(linedust[good].h_alpha[0],linedust[good].h_alpha[1],$
            hb[good].h_beta_dust_predict[0],hb[good].h_beta_dust_predict[1],/quotient)
          
          hb[good].ebv_hahb_predict = get_ebv(hahb,decrement_err=hahb_err,ebv_err=ebv_err) > 0.0
          hb[good].ebv_hahb_err_predict = ebv_err
          
       endif

       if (tag_exist(linenodust,'H_BETA_DUST_PREDICT') eq 0L) then $
         linenodust = struct_addtags(temporary(linenodust),hb)

       if (tag_exist(linedust,'H_BETA_DUST_PREDICT') eq 0L) then $
         linedust = struct_addtags(temporary(linedust),hb) ; also update the input structure

; H-GAMMA

       kl = k_lambda(4340.464,_extra=extra)
       
       hg = {$
         H_GAMMA_DUST_PREDICT:   [0.0,-2.0], $
         H_GAMMA_NODUST_PREDICT: [0.0,-2.0], $
         EBV_HAHG_PREDICT:        0.0, $
         EBV_HAHG_ERR_PREDICT:    0.0, $
         EBV_HBHG_PREDICT:        0.0, $
         EBV_HBHG_ERR_PREDICT:   -1.0}
         
       hg = replicate(hg,nspec)

       if (ngood ne 0L) then begin

          ebv = linenodust[good].ebv_hahb
          ebv_err = linenodust[good].ebv_hahb_err

          flux = linenodust[good].h_alpha[0] / R_HaHg
          ferr = linenodust[good].h_alpha[1] / R_HaHg
          
          hg[good].h_gamma_nodust_predict[0] = flux
          hg[good].h_gamma_nodust_predict[1] = ferr

          hg[good].h_gamma_dust_predict[0] = flux * 10^(-0.4*ebv*kl)
          hg[good].h_gamma_dust_predict[1] = sqrt( (ferr*10.0^(0.4*ebv*kl))^2.0 + $
            (flux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*ebv_err)^2.0 )

;         plot, alog10(hg[good].h_gamma_dust_predict[0]), linedust[good].h_gamma[0] / $
;           hg[good].h_gamma_dust_predict[0], xsty=3, ysty=3, ps=4
;         djs_oplot, !x.crange, [1,1], line=0, thick=2.0

          hahg = linedust[good].h_alpha[0]/hg[good].h_gamma_dust_predict[0]
          hahg_err = im_compute_error(linedust[good].h_alpha[0],linedust[good].h_alpha[1],$
            hg[good].h_gamma_dust_predict[0],hg[good].h_gamma_dust_predict[1],/quotient)
          
          hg[good].ebv_hahg_predict = get_ebv(hahg,decrement_err=hahg_err,ebv_err=ebv_err,/HaHg) > 0.0
          hg[good].ebv_hahg_err_predict = ebv_err
          
       endif

       good2 = where((linenodust.h_beta[1] gt 0.0) and (linenodust.ebv_hahb_err gt 0.0),ngood2)
       if (ngood2 ne 0L) then begin

          ebv = linenodust[good2].ebv_hahb
          ebv_err = linenodust[good2].ebv_hahb_err

          hbhg = linedust[good2].h_beta[0]/hg[good2].h_gamma_dust_predict[0]
          hbhg_err = im_compute_error(linedust[good2].h_beta[0],linedust[good2].h_beta[1],$
            hg[good2].h_gamma_dust_predict[0],hg[good2].h_gamma_dust_predict[1],/quotient)
          
          hg[good2].ebv_hbhg_predict = get_ebv(hbhg,decrement_err=hbhg_err,ebv_err=ebv_err,/HbHg) > 0.0
          hg[good2].ebv_hbhg_err_predict = ebv_err
          
       endif

       if (tag_exist(linenodust,'H_GAMMA_DUST_PREDICT') eq 0L) then $
         linenodust = struct_addtags(temporary(linenodust),hg)

       if (tag_exist(linedust,'H_GAMMA_DUST_PREDICT') eq 0L) then $
         linedust = struct_addtags(temporary(linedust),hg) ; also update the input structure

; H-DELTA       
       
       kl = k_lambda(4101.734,_extra=extra)
       
       hd = {H_DELTA_DUST_PREDICT: [0.0,-2.0], H_DELTA_NODUST_PREDICT: [0.0,-2.0]}
       hd = replicate(hd,nspec)

       if (ngood ne 0L) then begin

          ebv = linenodust[good].ebv_hahb
          ebv_err = linenodust[good].ebv_hahb_err

          flux = linenodust[good].h_alpha[0] / R_HaHd
          ferr = linenodust[good].h_alpha[1] / R_HaHd
          
          hd[good].h_delta_nodust_predict[0] = flux
          hd[good].h_delta_nodust_predict[1] = ferr

          hd[good].h_delta_dust_predict[0] = flux * 10^(-0.4*ebv*kl)
          hd[good].h_delta_dust_predict[1] = sqrt( (ferr*10.0^(0.4*ebv*kl))^2.0 + $
            (flux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*ebv_err)^2.0 )

;         plot, alog10(hd[good].h_delta_dust_predict[0]), linedust[good].h_delta[0] / $
;           hd[good].h_delta_dust_predict[0], xsty=3, ysty=3, ps=4
;         djs_oplot, !x.crange, [1,1], line=0, thick=2.0

       endif

       if (tag_exist(linenodust,'H_DELTA_DUST_PREDICT') eq 0L) and $
         (tag_exist(linenodust,'H_DELTA_NODUST_PREDICT') eq 0L) then $
         linenodust = struct_addtags(temporary(linenodust),hd)

       if (tag_exist(linedust,'H_DELTA_DUST_PREDICT') eq 0L) and $
         (tag_exist(linedust,'H_DELTA_NODUST_PREDICT') eq 0L) then $
         linedust = struct_addtags(temporary(linedust),hd) ; also update the input structure

    endif
       
; finally correct the LICK_HD_A index for contamination from H-delta
; emission  

    if tag_exist(linedust,'LICK_HD_A') and tag_exist(linedust,'H_DELTA') then begin

       good = where((linenodust.h_delta_continuum[1] gt 0.0) and $
         (linenodust.h_delta_dust_predict[1] gt 0.0),ngood)

       hdcor = replicate({LICK_HD_A_COR: [0.0,-2.0]},nspec)
       if (ngood ne 0L) then begin

          hd_dust_predict_EW = linenodust[good].h_delta_dust_predict[0] / linenodust[good].h_delta_continuum[0]
          hd_dust_predict_EW_err = im_compute_error(linenodust[good].h_delta_dust_predict[0],$
            linenodust[good].h_delta_dust_predict[1],linenodust[good].h_delta_continuum[0],$
            linenodust[good].h_delta_continuum[1],/quotient)

; *add* the correction to make the absorption indices bigger
          
          hdcor[good].lick_hd_a_cor[0] = linedust[good].lick_hd_a[0] + hd_dust_predict_EW
          hdcor[good].lick_hd_a_cor[1] = sqrt(linedust[good].lick_hd_a[1]^2 + hd_dust_predict_EW_err^2)

       endif

       if (tag_exist(linenodust,'LICK_HD_A_COR') eq 0L) then $
         linenodust = struct_addtags(temporary(linenodust),hdcor)

       if (tag_exist(linedust,'LICK_HD_A_COR') eq 0L) then $
         linedust = struct_addtags(temporary(linedust),hdcor) ; also update the input structure
          
    endif
    
return, linenodust
end    
