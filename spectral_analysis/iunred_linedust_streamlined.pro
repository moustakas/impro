;+
; NAME:
;       IUNRED_LINEDUST_STREAMLINED()
;
; PURPOSE:
;       Correct emission-line fluxes for dust extinction using the
;       observed Balmer decrement.
;
; INPUTS:
;       linedust - input data structure
;
; OPTIONAL INPUTS:
;       snrcut_linedust - 
;       cutsig          - 
;       extra           - keywords for GET_EBV()
;
; KEYWORD PARAMETERS:
;       HaHb        - de-redden using the H-alpha/H-beta line ratio 
;                     (default) 
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
;       J. Moustakas, 2008 Feb 20, NYU - excised from IUNRED_LINEDUST
;                                        and made simple!
;
; Copyright (C) 2008, John Moustakas
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

function iunred_linedust_streamlined, linedust, snrcut_linedust=snrcut_linedust, $
  cutsig=cutsig, HaHb=HaHb, nopropagate=nopropagate, allow_r_min=allow_r_min, $
  silent=silent, _extra=extra
    
    nspec = n_elements(linedust)
    if (nspec eq 0L) then begin
       doc_library, 'iunred_linedust_streamlined'
       return, -1L
    endif

    if (n_elements(snrcut_linedust) eq 0L) then snrcut_linedust = 5.0
    if (n_elements(cutsig) eq 0L) then cutsig = 1.0

    if not keyword_set(silent) then splog, 'S/N > '+string(snrcut_linedust,format='(G0.0)')+'.'

    linename = strlowcase(strcompress(linedust[0].linename,/remove))
    nline = n_elements(linename)

    R_HaHb = 2.86     ; intrinsic Ha/Hb ratio
    R_HaHb_min = 2.75 ; minimum physically plausible Balmer decrements
    
; initialize and fill the output structure

    newtags = {$
      hahb:  0.0, hahb_err:  -1.0, ebv_hahb: 0.0, ebv_hahb_err: -1.0, ebv_hahb_negative: 0B, $
      ehbha: 0.0, ehbha_err: -1.0}
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

; reddening values    
    
    HaHb = 1L
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

; the T structure data points can be corrected for reddening but the
; NORED structure data point errors must be set to -1.0
    
    t = linenodust[good]
    if (nbad ne 0L) then nored = linenodust[bad]

; loop on each line and correct for dust.  also, if emission line
; luminosities have been computed then correct those for dust as well;
; also correct the emission-line continuum and form new
; "dust-corrected" equivalent widths (in the general case where RATIO
; is not equal to 1.0)

    tags = tag_names(linenodust)
    for i = 0L, nline-1L do begin

; only de-redden well-measured emission lines; should the upper limits
; be de-reddened as well?  don't for now.

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
       t.(match) = float(nodust)
       if (nbad ne 0L) then nored.(match) = [0.0,-1.0]

    endfor
       
; now add the line flux corrections to LINENODUST

    linenodust[good] = struct_trimtags(t,select_tags=tag_names(t))
    if (nbad ne 0L) then linenodust[bad] = struct_trimtags(nored,select_tags=tag_names(nored))

    
return, linenodust
end    
