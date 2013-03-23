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
;       snrcut_linedust - minimum S/N on H-alpha and H-beta (default 5)
;       use_HaHb        - intrinsic Ha/Hb ratio to use (default 2.86)
;       cutsig          - number of standard deviations that the
;                         observed Ha/Hb ratio can deviate from the
;                         minimum physical intrinsic ratio (default 1;
;                         see COMMENTS)
;       extra           - keywords specifying the attenuation curve in
;                         GET_EBV() (default O'Donnell 1994)
;
; KEYWORD PARAMETERS:
;       nopropagate - do not propagate the error in the reddening
;                     correction to the unreddened fluxes (not
;                     recommended!) 
;       allow_r_min - use the minimum physically plausible Balmer
;                     decrement, rather than the mean (T=10^4 K) value 
;       silent      - suppress messages to STDOUT
;
; OUTPUTS:
;       linenodust - output data structure in which all the
;                    line-fluxes have been reddening-corrected, and
;                    the reddening information is stored in additional
;                    tags 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Upper limits on the Balmer lines are not treated. 
;
;       Do not compute E(B-V) for objects with unphysically low Balmer
;       decrements; they must be within CUTSIG sigma of the
;       adopted/minimum decrement; there is no upper limit on the
;       Balmer decrement (dusty!); note, however, that objects can
;       have negative E(Hb-Ha) that is zero within the errors;
;       however, for these objects keep the measured Balmer decrements
;       (ie, in general they will be less than 2.86); jm04jul22uofa:
;       for objects with E(B-V)<0 set E(B-V)=0 and flag them;
;       jm05feb03uofa: new conservative policy: require the decrement
;       to be within CUTSIG of the **adopted** intrinsic decrement,
;       **not** the minimum possible decrement; jm08feb06nyu: use a
;       keyword to allow the smallest possible decrement 
;
;       Note that the code *could* take into account the error in the
;       intrinsic HaHb ratio (due to, e.g., ignorance of the
;       appropriate electron temperature/density), but for now we
;       assume the uncertainty is negligible.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2008 Mar 07, NYU - simplified considerably and
;          made much, much faster; errors now computed correctly,
;          accounting for the covariance between the error in the
;          Balmer decrement and the reddening-corrected H-alpha and
;          H-beta line-fluxes
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

function iunred_linedust, linedust, snrcut_linedust=snrcut_linedust, $
  use_HaHb=use_HaHb, cutsig=cutsig, nopropagate=nopropagate, $
  allow_r_min=allow_r_min, silent=silent, _extra=extra
    
    nspec = n_elements(linedust)
    if (nspec eq 0L) then begin
       doc_library, 'iunred_linedust'
       return, -1L
    endif

    if (n_elements(snrcut_linedust) eq 0L) then snrcut_linedust = 5.0
    if (n_elements(cutsig) eq 0L) then cutsig = 1.0

    if not keyword_set(silent) then splog, 'S/N > '+$
      string(snrcut_linedust,format='(G0.0)')

; verify that the input data structure is compatible

    if (tag_exist(linedust[0],'LINENAME') eq 0L) then begin
       splog, 'LINEDUST is missing required LINENAME tag'
       return, linedust
    endif
    
    linename = strlowcase(strcompress(linedust[0].linename,/remove))
    nline = n_elements(linename)

; initialize and fill the output structure; R_HaHb_min is the minimum
; physically plausible Balmer decrement 

    R_HaHb_err = 0.0
    R_HaHb_min = 2.75 

    if (n_elements(use_HaHb) ne 0L) then R_HaHb = use_HaHb else R_HaHb = 2.86
    if keyword_set(allow_r_min) then this_R_HaHb = R_HaHb_min else this_R_HaHb = R_HaHb

    newtags = {hahb:  0.0, hahb_err:  -1.0, ebv_hahb: 0.0, $
      ebv_hahb_err: -1.0, ehbha: 0.0, ehbha_err: -1.0}
    newtags = replicate(newtags,nspec)

    linenodust = struct_addtags(linedust,temporary(newtags))

; ---------------------------------------------------------------------------
; Ha/Hb ratio, E(Hb-Ha) value, and corresponding E(B-V)
; ---------------------------------------------------------------------------

    if tag_exist(linedust,'H_ALPHA') and tag_exist(linedust,'H_BETA') then begin
       
; keep all objects with well-measured H-alpha and H-beta lines

       ha = linenodust.h_alpha[0]
       haerr = linenodust.h_alpha[1]
       hb = linenodust.h_beta[0]
       hberr = linenodust.h_beta[1]
       
       goodsnr = where((haerr gt 0.0) and (ha/haerr gt snrcut_linedust) and $
         (hberr gt 0.0) and (hb/hberr gt snrcut_linedust),ngoodsnr)

       if (ngoodsnr ne 0L) then begin

          linenodust[goodsnr].hahb = ha[goodsnr]/hb[goodsnr]
          linenodust[goodsnr].hahb_err = sqrt((haerr[goodsnr]/hb[goodsnr])^2.0 + $
            hberr[goodsnr]^2.0*(ha[goodsnr]/hb[goodsnr]^2.0)^2.0)
          
          keep = where((linenodust[goodsnr].hahb+cutsig*$ ; account for the error in HaHb
            linenodust[goodsnr].hahb_err ge this_R_HaHb),nkeep)

          if (nkeep ne 0L) then begin

             linenodust[goodsnr[keep]].ebv_hahb = get_ebv(linenodust[goodsnr[keep]].hahb,$
               decrement_err=linenodust[goodsnr[keep]].hahb_err,ebv_err=ebv_err,color=color,$
               err_color=color_err,/HaHb,_extra=extra,use_HaHb=use_HaHb)

             linenodust[goodsnr[keep]].ebv_hahb_err = ebv_err
             linenodust[goodsnr[keep]].ehbha = color
             linenodust[goodsnr[keep]].ehbha_err = color_err

          endif

          neg = where((linenodust.ebv_hahb lt 0.0),nneg)
          good = where((linenodust.ebv_hahb_err gt 0.0),ngood)

          if (nneg ne 0L) then begin
             if (not keyword_set(silent)) then splog, 'Identified '+$
               string(nneg,format='(I0)')+'/'+string(ngood,format='(I0)')+$
               ' [of '+string(nspec,format='(I0)')+$
               ' objects] with negative E(B-V) [Ha/Hb].'
             linenodust[neg].ebv_hahb = 0.0
          endif 
          
       endif else splog, 'WARNING: No good Ha/Hb Balmer decrement measurements!'

    endif

; loop on each line and correct for dust
    tags = tag_names(linenodust)
    hamatch = where('h_alpha' eq strlowcase(tags),nhamatch)
    hbmatch = where('h_beta' eq strlowcase(tags),nhbmatch)
    
    for iline = 0, nline-1 do begin

       match = where(linename[iline] eq strlowcase(tags),nmatch)
       if (nmatch eq 0L) then begin
          splog, 'Emission line '+strupcase(linename[iline])+' not found...continuing.'
          continue
       endif

       wavematch = where(linename[iline]+'_wave' eq $
         strlowcase(tags),nwavematch)
       if (nwavematch eq 0L) then begin
          splog, 'Emission line wavelength '+$
            strupcase(linename[iline])+'_WAVE not found...continuing'
          continue
       endif

       lineflux = reform(linenodust.(match)[0,*])
       lineferr = reform(linenodust.(match)[1,*])
       linewave = linenodust.(wavematch) ; [Angstrom]

       newflux = lineflux*0.0
       newferr = lineferr*0.0-1.0

       nice = where((linenodust.ebv_hahb_err gt 0.0) and (lineferr gt 0.0),nnice)

       if (nnice ne 0L) then begin
; correct for dust
          ebv = linenodust[nice].ebv_hahb
          ebv_err = linenodust[nice].ebv_hahb_err
          kl = k_lambda(linewave[nice],_extra=extra)
          newflux[nice] = lineflux[nice]*10.0^(0.4*ebv*kl)
; now compute the uncertainty
          if keyword_set(nopropagate) then begin
             newferr[nice] = lineferr[nice]*10.0^(0.4*ebv*kl)
          endif else begin
; take into account the covariance between the observed Balmer
; decrement, and the emission line under consideration
             mlam = -kl/(k_lambda(4861.33,_extra=extra)-k_lambda(6562.80,_extra=extra))
             ha = reform((linedust[nice].(hamatch))[0,*])
             hb = reform((linedust[nice].(hbmatch))[0,*])
             haerr = reform((linedust[nice].(hamatch))[1,*])
             hberr = reform((linedust[nice].(hbmatch))[1,*])
; case 1a - this line is neither H-alpha nor H-beta; compute the
; "relative variance"
             if (match ne hamatch) and (match ne hbmatch) then begin 
                newrelvar = (lineferr[nice]/lineflux[nice])^2.0 + $
                  mlam^2.0*((haerr/ha)^2.0 + (hberr/hb)^2.0 + (R_HaHb_err/R_HaHb)^2.0)
             endif
; case 1b - this line is H-alpha
             if (match eq hamatch) then begin 
                newrelvar = (1-mlam)^2.0*(haerr/ha)^2.0 + mlam^2.0*((hberr/hb)^2.0 + $
                  (R_HaHb_err/R_HaHb)^2.0)
             endif
; case 1c - this line is H-beta
             if (match eq hbmatch) then begin 
                newrelvar = (1-mlam)^2.0*(hberr/hb)^2.0 + mlam^2.0*((haerr/ha)^2.0 + $
                  (R_HaHb_err/R_HaHb)^2.0)
             endif
             newferr[nice] = sqrt(newflux[nice]^2.0*newrelvar)
; old bit of code
;            tempflux = dust_correct(lineflux[nice],linewave[nice],$
;              err_lineflux=lineferr[nice],ebv=ebv,err_ebv=ebv_err,$
;              err_newflux=tempferr)
;;           newflux[nice] = tempflux & newferr[nice] = tempferr
          endelse

       endif 

       linenodust.(match) = transpose([ [newflux], [newferr] ])
       
    endfor

return, linenodust
end    
