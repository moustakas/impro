;+
; NAME:
;       IM_ASSIGN_R23BRANCH()
;
; PURPOSE:
;       Given the output from IM_ABUNDANCE() assign the appropriate
;       R23 branch. 
;
; INPUTS:
;       abund - output structure from IM_ABUNDANCE() or MZ_ABUNDANCE() 
;
; OPTIONAL INPUTS:
;       mb           - absolute B-band magnitude (if BRANCHMETHOD=1)
;       branchmethod - technique to use to assign branches; in either case,
;                      we average the lower- and upper-branch
;                      abundances of ambigious objects that are in the
;                      turn-around region  
;          0L - [NII]/Ha<-1 --> lower, [NII]/Ha>-1 --> upper
;          1L - use the luminosity-metallicity relation iteratively
;               (requires MB input); only recommended if most objects
;               are on the upper or lower branch
;          2L - use 
;          3L - use [N II]/Ha and [N II]/[O II], following Contini et
;               al. (2002)
;          4L - put them all on the upper branch!
;          5L - assign them according to R23BRANCH
;       subindx - only assign branches to ABUND[SUBINDX] (see, e.g.,
;                 WRITE_AGES_MZ_SAMPLE)
;
; KEYWORD PARAMETERS:
;       kk04 - consider the Kobulnicky & Kewley (2004) abundances
;              (default) 
;       pt05 - consider the Pilyugin & Thuan (2005) abundances
;       m91  - consider the McGaugh (1991) abundances
;
; OUTPUTS:
;       result - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 18, U of A - written
;       jm07sep18nyu - added JUSTEW and JUSTFLUX keywords
;       jm08feb11nyu - added FRACCUT optional input and SILENT keyword
;       jm08oct22nyu - FRACCUT is now obsolete because the formal
;         1-sigma calculation is done (see SINGS paper)
;
; Copyright (C) 2006-2008, John Moustakas
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

function im_assign_r23branch, abund, mb=mb, branchmethod=branchmethod, $
  subindx=subindx, r23branch=r23branch1, o32cut=o32cut, n2cut=n2cut, $
  fraccut=fraccut, logr23cut=logr23cut, justew=justew, justflux=justflux, kk04=kk04, $
  pt05=pt05, m91=m91, test=test, debug=debug, silent=silent

    nobj = n_elements(abund)
    if (nobj eq 0L) then begin
       doc_library, 'im_assign_r23branch'
       return, -1L
    endif
    
    if (n_elements(subindx) eq 0L) then subindx = lindgen(nobj)
    if (n_elements(o32cut) eq 0L) then o32cut = -0.2
    if (n_elements(n2cut) eq 0L) then n2cut = -1.1
    if (n_elements(fraccut) eq 0L) then fraccut = 0.2
    if (n_elements(logr23cut) eq 0L) then logr23cut = 0.8

    if (n_elements(branchmethod) eq 0L) then branchmethod = 0L
    if (branchmethod eq 1L) then begin
       if (n_elements(mb) eq 0L) then begin
          splog, 'MB input required for BRANCHMETHOD=1!'
          return, -1L
       endif
       if (n_elements(mb) ne nobj) then begin
          splog, 'ABUND and MB have incompatible dimensions.'
          return, -1L
       endif
    endif

    if (branchmethod eq 5L) then begin
       if (n_elements(r23branch1) eq 0L) then begin
          splog, 'R23BRANCH input required for BRANCHMETHOD=5!'
          return, -1L
       endif else r23branch = r23branch1
       if (n_elements(r23branch1) ne nobj) then begin
          if (n_elements(r23branch1) eq 1L) then r23branch = replicate(r23branch1,nobj) else begin
             splog, 'ABUND and R23BRANCH have incompatible dimensions.'
             return, -1L
          endelse
       endif
    endif

    if (n_elements(kk04) eq 0L) and (n_elements(pt05) eq 0L) and (n_elements(m91) eq 0L) then kk04 = 1L
    
    if keyword_set(kk04) then begin
       suffix = 'KK04' & newsuffix = suffix
    endif
    if keyword_set(pt05) then begin
       suffix = 'PT05' & newsuffix = suffix
    endif
    if keyword_set(m91) then begin
       suffix = 'M91' & newsuffix = suffix
    endif

; default: compute abundances based on both EWS and FLUXES 
    
    if (not keyword_set(justew)) and (not keyword_set(justflux)) then begin
       justew = 1L
       justflux = 1L
    endif
    
    allresult = {$
      ewalpha:                1.0,$
      r23branch:              '?',$
      r23branch_ew:           '?',$
      zstrong_niiha:       -999.0,$
      zstrong_niiha_err:   -999.0,$
      zstrong_niioii:      -999.0,$
      zstrong_niioii_err:  -999.0,$
      zstrong_o32:         -999.0,$
      zstrong_o32_err:     -999.0,$
      zstrong_ew_o32:      -999.0,$
      zstrong_ew_o32_err:  -999.0,$
      zstrong_r23:         -999.0,$
      zstrong_r23_err:     -999.0,$
      zstrong_ew_r23:      -999.0,$
      zstrong_ew_r23_err:  -999.0,$
      zstrong_p:           -999.0,$
      zstrong_p_err:       -999.0,$
      zstrong_ew_p:        -999.0,$
      zstrong_ew_p_err:    -999.0,$
      zstrong_logu:        -999.0,$
      zstrong_logu_err:    -999.0,$
      zstrong_ew_logu:     -999.0,$
      zstrong_ew_logu_err: -999.0,$
      zstrong_12oh:        -999.0,$
      zstrong_12oh_err:    -999.0,$
      zstrong_ew_12oh:     -999.0,$
      zstrong_ew_12oh_err: -999.0}
    allresult = replicate(allresult,nobj)

    struct_assign, abund, allresult, /nozero ; this has to go before the next line
    result = allresult[subindx]              ; select a subset of objects
    ngalaxy = n_elements(result)

; assign the appropriate suffix to the output tags    

    finaltags = tag_names(allresult[0])
    alter = where(strmatch(finaltags,'*12oh*',/fold) or strmatch(finaltags,'*logu*',/fold) or strmatch(finaltags,'*r23branch*',/fold))
    finaltags[alter] = repstr(finaltags[alter]+'_'+newsuffix,'_ERR_'+newsuffix,'_'+newsuffix+'_ERR')

;   finaltags = ['r23branch_'+newsuffix+['','_ew'],'zstrong_12oh_'+$ ; final output tags
;     newsuffix+['','_err'],'zstrong_ew_12oh_'+newsuffix+['','_err']]

    if (n_elements(mb) ne 0L) then mbsub = mb[subindx]
    if (n_elements(r23branch) ne 0L) then r23branchsub = r23branch[subindx]
    nsubindx = n_elements(subindx)
    
    data = struct_trimtags(abund[subindx],select=['ZSTRONG_12OH_'+suffix+'_*','ZSTRONG_EW_12OH_'+suffix+'_*',$
      'ZSTRONG_LOGU_'+suffix+'_*','ZSTRONG_EW_LOGU_'+suffix+'_*'])
    data = im_struct_trimtags(data,select=tag_names(data),$
      newtags=repstr(tag_names(data),suffix+'_',''))

    tremonti_coeff = [5.276,-0.186]
    mbaxis = findgen(((-10.0)-(-30.0))/0.01)*0.01+(-30.0)

    maxiter = 5L
    plotsym, 0, 0.4, /fill
    
; #########################
; fluxes    
; #########################

    if keyword_set(justflux) then begin
       
       if (not keyword_set(silent)) then splog, 'Assigning R23 branches from fluxes:'

       case branchmethod of

          0L: begin             ; just use [N II]/Ha
             
             lo = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha le n2cut) and $
               (data.zstrong_12oh_lower gt -900) and (data.zstrong_12oh_upper gt -900) and $
               (data.zstrong_12oh_lower lt data.zstrong_12oh_upper),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch        = 'L_N2'
                result[lo].zstrong_12oh     = data[lo].zstrong_12oh_lower
                result[lo].zstrong_12oh_err = data[lo].zstrong_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_logu     = data[lo].zstrong_logu_lower
                   result[lo].zstrong_logu_err = data[lo].zstrong_logu_lower_err
                endif
             endif
             
             up = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha gt n2cut) and $
               (data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
               (data.zstrong_12oh_upper gt data.zstrong_12oh_lower),nup)
             if (nup ne 0L) then begin
                result[up].r23branch        = 'U_N2'
                result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
                result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_logu     = data[up].zstrong_logu_upper
                   result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
                endif
             endif

          end

          1L: begin             ; iterative LZ method
             
; begin by assigning branches using the O32 ratio (or put everything
; on the upper branch); on subsequent iterations, use the abundance
; that is closest to the abundance predicted by the LZ relation
; (either derived from the data, or based on Tremonti et al. 2004)

             good = where((data.zstrong_12oh_upper gt -900.0) and (data.zstrong_o32 gt -900.0) and (mbsub gt -900.0),ngood)
             if (ngood ne 0L) then begin

                mbgood = mbsub[good]
                ohgood = mbgood*0.0
                
                up = where((data[good].zstrong_o32 le o32cut),nup,comp=lo,ncomp=nlo)
                if (nup ne 0L) then ohgood[up] = data[good[up]].zstrong_12oh_upper
                if (nlo ne 0L) then ohgood[lo] = data[good[lo]].zstrong_12oh_lower
;               ohgood = data[good].zstrong_12oh_upper
                ohgood_err = ohgood*0.0

                for ii = 0L, maxiter-1L do begin
                   
                   lzoh = poly(mbgood,tremonti_coeff)
;                  sixlin, mbgood, ohgood, lzint, siga, lzslope, sigb
;                  lzoh = poly(mbgood,[lzint[2],lzslope[2]])
;                  if keyword_set(debug) then begin
;                     plot, mbgood, ohgood, xsty=3, ysty=3, ps=3, $
;                       title='Flux/'+strupcase(suffix), charsize=1.5, charthick=2.0, $
;                       xtitle=textoidl('M_{B}'), ytitle='12+log(O/H)'
;                     oplot, mbgood, lzoh, line=0, thick=2
;                     cc = get_kbrd(1)
;                  endif

                   up = where((abs(lzoh-data[good].zstrong_12oh_upper) lt abs(lzoh-data[good].zstrong_12oh_lower)),$
                     nup,comp=lo,ncomp=nlo)

                   if (nup ne 0L) then begin
                      result[good[up]].r23branch        = 'U_LZ'
                      result[good[up]].zstrong_12oh     = data[good[up]].zstrong_12oh_upper
                      result[good[up]].zstrong_12oh_err = data[good[up]].zstrong_12oh_upper_err
                      if keyword_set(kk04) then begin
                         result[good[up]].zstrong_logu     = data[good[up]].zstrong_logu_upper
                         result[good[up]].zstrong_logu_err = data[good[up]].zstrong_logu_upper_err
                      endif
                   endif
                   if (nlo ne 0L) then begin
                      result[good[lo]].r23branch        = 'L_LZ'
                      result[good[lo]].zstrong_12oh     = data[good[lo]].zstrong_12oh_lower
                      result[good[lo]].zstrong_12oh_err = data[good[lo]].zstrong_12oh_lower_err
                      if keyword_set(kk04) then begin
                         result[good[lo]].zstrong_logu     = data[good[lo]].zstrong_logu_lower
                         result[good[lo]].zstrong_logu_err = data[good[lo]].zstrong_logu_lower_err
                      endif
                   endif

                   ohgood = result[good].zstrong_12oh ; update the abundances
                   ohgood_err = result[good].zstrong_12oh_err

                endfor
                
             endif
             
          end

          2L: begin             ; combination of [N II]/Ha and LZ assignment

; assign the branch according to [N II]/Ha
             
             lo = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha le n2cut) and $
               (data.zstrong_12oh_lower gt -900) and (data.zstrong_12oh_upper gt -900) and $
               (data.zstrong_12oh_lower lt data.zstrong_12oh_upper),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch        = 'L_N2'
                result[lo].zstrong_12oh     = data[lo].zstrong_12oh_lower
                result[lo].zstrong_12oh_err = data[lo].zstrong_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_logu     = data[lo].zstrong_logu_lower
                   result[lo].zstrong_logu_err = data[lo].zstrong_logu_lower_err
                endif
             endif
             
             up = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha gt n2cut) and $
               (data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
               (data.zstrong_12oh_upper gt data.zstrong_12oh_lower),nup)
             if (nup ne 0L) then begin
                result[up].r23branch        = 'U_N2'
                result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
                result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_logu     = data[up].zstrong_logu_upper
                   result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
                endif
             endif

; fit a bisector to the data (or use the Tremonti LZ relation)          
             
             good = where((result.zstrong_12oh gt -900.0) and (mbsub gt -900.0),ngood)
             if (ngood ne 0L) then begin
;               sixlin, mbsub[good], result[good].zstrong_12oh, lzint, siga, lzslope, sigb
;               coeff = [lzint[2],lzslope[2]]
                coeff = tremonti_coeff
             endif 

; assign branches to the ambigious objects          

;            need = lindgen(nsubindx) & nneed = nsubindx
;            need = where((result.zstrong_12oh lt -900.0) and (data.zstrong_12oh_upper gt -900) and $
;              (data.zstrong_12oh_lower gt -900),nneed)
             need = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900),nneed)
             if (nneed ne 0L) then begin

                lzoh = poly(mbsub[need],coeff)
                up = where((abs(lzoh-data[need].zstrong_12oh_upper) lt abs(lzoh-data[need].zstrong_12oh_lower)),$
                  nup,comp=lo,ncomp=nlo)

                if (nup ne 0L) then begin
                   result[need[up]].r23branch        = 'U_LZ'
                   result[need[up]].zstrong_12oh     = data[need[up]].zstrong_12oh_upper
                   result[need[up]].zstrong_12oh_err = data[need[up]].zstrong_12oh_upper_err
                   if keyword_set(kk04) then begin
                      result[need[up]].zstrong_logu     = data[need[up]].zstrong_logu_upper
                      result[need[up]].zstrong_logu_err = data[need[up]].zstrong_logu_upper_err
                   endif
                endif
                if (nlo ne 0L) then begin
                   result[need[lo]].r23branch        = 'L_LZ'
                   result[need[lo]].zstrong_12oh     = data[need[lo]].zstrong_12oh_lower
                   result[need[lo]].zstrong_12oh_err = data[need[lo]].zstrong_12oh_lower_err
                   if keyword_set(kk04) then begin
                      result[need[lo]].zstrong_logu     = data[need[lo]].zstrong_logu_lower
                      result[need[lo]].zstrong_logu_err = data[need[lo]].zstrong_logu_lower_err
                   endif
                endif

             endif
             
          end 

          3L: begin             ; use [N II]/Ha and [N II]/[O II]
             
             lo = where((data.zstrong_niiha gt -900.0) and (data.zstrong_niiha lt -1.0) and $
               (data.zstrong_niioii gt -900.0) and (data.zstrong_niioii lt -1.05) and $
               (data.zstrong_12oh_lower gt -900.0) and (data.zstrong_12oh_upper gt -900.0) and $
               (data.zstrong_12oh_lower lt data.zstrong_12oh_upper),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch        = 'L_N2_N2O2'
                result[lo].zstrong_12oh     = data[lo].zstrong_12oh_lower
                result[lo].zstrong_12oh_err = data[lo].zstrong_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_logu     = data[lo].zstrong_logu_lower
                   result[lo].zstrong_logu_err = data[lo].zstrong_logu_lower_err
                endif
             endif
             
             up = where((data.zstrong_niiha gt -900.0) and (data.zstrong_niiha gt -1.0) and $
               (data.zstrong_niioii gt -900.0) and (data.zstrong_niioii gt -0.8) and $
               (data.zstrong_12oh_upper gt -900.0) and (data.zstrong_12oh_lower gt -900.0) and $
               (data.zstrong_12oh_upper gt data.zstrong_12oh_lower),nup)
             if (nup ne 0L) then begin
                result[up].r23branch        = 'U_N2_N2O2'
                result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
                result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_logu     = data[up].zstrong_logu_upper
                   result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
                endif
             endif

          end

          4L: begin             ; put them all on the upper branch

             good = where(data.zstrong_12oh_upper gt -900.0,ngood)
             if (ngood ne 0L) then begin
                result[good].r23branch        = 'U'
                result[good].zstrong_12oh     = data[good].zstrong_12oh_upper
                result[good].zstrong_12oh_err = data[good].zstrong_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[good].zstrong_logu     = data[good].zstrong_logu_upper
                   result[good].zstrong_logu_err = data[good].zstrong_logu_upper_err
                endif
             endif

          end

          5L: begin             ; assign them according to R23BRANCH

             up = where((data.zstrong_12oh_upper gt -900.0) and (r23branchsub eq 'U'),nup)
             if (nup ne 0L) then begin
                result[up].r23branch        = 'U'
                result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
                result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_logu     = data[up].zstrong_logu_upper
                   result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
                endif
             endif

             lo = where((data.zstrong_12oh_lower gt -900.0) and (r23branchsub eq 'L'),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch        = 'L'
                result[lo].zstrong_12oh     = data[lo].zstrong_12oh_lower
                result[lo].zstrong_12oh_err = data[lo].zstrong_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_logu     = data[lo].zstrong_logu_lower
                   result[lo].zstrong_logu_err = data[lo].zstrong_logu_lower_err
                endif
             endif

          end

          else: begin
             splog, 'BRANCHMETHOD value not supported!'
             return, -1L
          endelse

       endcase

; adopt the "average" lower/upper abundance of objects with ambiguous
; abundances, near the turn-around region; reject objects with a
; "frac" value less than fraccut (see IM_ABUNDANCE for details)

;      ambig = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
;        ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower) and $ ; off the calibration
;        (data.zstrong_12oh_upper+data.zstrong_12oh_upper_err) lt $
;        (data.zstrong_12oh_lower-data.zstrong_12oh_lower_err)) or $
;        ((data.zstrong_12oh_upper gt data.zstrong_12oh_lower) and $ ; on the R23 calibration
;        (data.zstrong_12oh_upper-data.zstrong_12oh_upper_err) lt $
;        (data.zstrong_12oh_lower+data.zstrong_12oh_lower_err)),nambig)
       ambig = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
;        (abund.zstrong_r23 gt logr23cut) and $
         ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower)) or $ ; off the R23 calibration
         ((data.zstrong_12oh_upper gt data.zstrong_12oh_lower) and $ ; on the R23 calibration, but ambiguous
         (data.zstrong_12oh_upper-data.zstrong_12oh_upper_err) lt $
         (data.zstrong_12oh_lower+data.zstrong_12oh_lower_err)),nambig)
;      ambig = where((data.zstrong_12oh_upper gt -900.0) and (data.zstrong_12oh_lower gt -900.0) and $
;        (data.zstrong_12oh_upper lt data.zstrong_12oh_lower),nambig)
       if (not keyword_set(silent)) then splog, '   Ambiguous: (O/H)_lower>(O/H)_upper: '+$
         string(nambig,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngalaxy),format='(F12.1)'),2)+'%).'

; reject objects whose upper- and lower-branch abundances differ by
; more than 1-sigma, where 1-sigma is given by the width of the
; overlapping histogram (12oh_avg)
       if (nambig ne 0L) then begin
          good = where((data[ambig].zstrong_12oh_avg+data[ambig].zstrong_12oh_avg_err gt data[ambig].zstrong_12oh_lower) and $
            (data[ambig].zstrong_12oh_avg-data[ambig].zstrong_12oh_avg_err lt data[ambig].zstrong_12oh_upper),ngood,$
            comp=reject,ncomp=nreject)
          if (not keyword_set(silent)) then splog, '     Retain: '+$
            string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%).'
          if (ngood ne 0L) then begin
             result[ambig[good]].r23branch        = 'A'
             result[ambig[good]].zstrong_12oh     = data[ambig[good]].zstrong_12oh_avg
             result[ambig[good]].zstrong_12oh_err = data[ambig[good]].zstrong_12oh_avg_err
             if keyword_set(kk04) then begin
                result[ambig[good]].zstrong_logu     = data[ambig[good]].zstrong_logu_avg
                result[ambig[good]].zstrong_logu_err = data[ambig[good]].zstrong_logu_avg_err
             endif
          endif
          if (not keyword_set(silent)) then splog, '     Reject: '+$
            string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*nreject/float(nambig),format='(I0)'),2)+'%).'
          if (nreject ne 0L) then begin
             result[ambig[reject]].r23branch        = 'Rejected'
             result[ambig[reject]].zstrong_12oh     = -999.0
             result[ambig[reject]].zstrong_12oh_err = -999.0
             result[ambig[reject]].zstrong_logu     = -999.0
             result[ambig[reject]].zstrong_logu_err = -999.0
          endif
       endif
; OLD CODE THAT USES FRACCUT
;      if (nambig ne 0L) then begin
;         good = where((data[ambig].zstrong_12oh_frac gt fraccut),ngood,comp=reject,ncomp=nreject)
;         if (ngood ne 0L) then begin
;            if (not keyword_set(silent)) then splog, '     Retain (F>'+$
;              strtrim(string(fraccut,format='(F12.1)'),2)+'): '+$
;              string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
;              strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%).'
;            result[ambig[good]].r23branch        = 'A'
;            result[ambig[good]].zstrong_12oh     = data[ambig[good]].zstrong_12oh_avg
;            result[ambig[good]].zstrong_12oh_err = data[ambig[good]].zstrong_12oh_avg_err
;            if keyword_set(kk04) then begin
;               result[ambig[good]].zstrong_logu     = data[ambig[good]].zstrong_logu_avg
;               result[ambig[good]].zstrong_logu_err = data[ambig[good]].zstrong_logu_avg_err
;            endif
;         endif
;         if (nreject ne 0L) then begin
;            if (not keyword_set(silent)) then splog, '     Reject (F<'+$
;              strtrim(string(fraccut,format='(F12.1)'),2)+'): '+$
;              string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
;              strtrim(string(100.0*nreject/float(nambig),format='(I0)'),2)+'%).'
;            result[ambig[reject]].r23branch        = 'Rejected'
;            result[ambig[reject]].zstrong_12oh     = -999.0
;            result[ambig[reject]].zstrong_12oh_err = -999.0
;            result[ambig[reject]].zstrong_logu     = -999.0
;            result[ambig[reject]].zstrong_logu_err = -999.0
;         endif
;      endif

; MORE OLD CODE       
;      nambig = 0L
;      if (nambig ne 0L) then begin
;         result[ambig].r23branch = 'A'
;         for i = 0L, nambig-1L do begin
;            oh     = [data[ambig[i]].zstrong_12oh_upper,data[ambig[i]].zstrong_12oh_lower]
;            oh_err = [data[ambig[i]].zstrong_12oh_upper_err,data[ambig[i]].zstrong_12oh_lower_err]
;            result[ambig[i]].zstrong_12oh     = total(oh/oh_err^2)/total(1.0/oh_err^2)
;            result[ambig[i]].zstrong_12oh_err = 1.0/sqrt(total(1.0/oh_err^2))
;         endfor
;      endif

; reject objects having O/H_upper<O/H_lower that are not within
; one-sigma of one another (these usually have R23>10) (this step
; should not be necessary)
       
;      reject = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
;        ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower) and $ 
;         ((data.zstrong_12oh_upper+data.zstrong_12oh_upper_err) lt $ ; note!
;          (data.zstrong_12oh_lower-data.zstrong_12oh_lower_err))),nreject)
;      if (nreject ne 0L) then begin
;         splog, 'Rejecting (flux): (O/H)_lower>(O/H)_upper: '+string(nreject,format='(I0)')+'/'+$
;           string(ngalaxy,format='(I0)')+' ('+string(round(100.0*nreject/float(ngalaxy)),format='(I0)')+'%).'
;         result[reject].r23branch        = 'Rejected'
;         result[reject].zstrong_12oh     = -999.0
;         result[reject].zstrong_12oh_err = -999.0
;      endif

       if keyword_set(debug) and (n_elements(mb) ne 0L) then begin
          good = where((mbsub gt -900.0) and (result.zstrong_12oh gt -900.0))
          plot, mbsub[good], result[good].zstrong_12oh, ps=8, xsty=3, ysty=3, $
            title='Flux/'+strupcase(suffix)+' - Final', charsize=1.5, charthick=2.0, $
            xtitle=textoidl('M_{B}'), ytitle='12+log(O/H)'
          up = where(strmatch(result[good].r23branch,'*U*') eq 1B,nup) 
          lo = where(strmatch(result[good].r23branch,'*L*') eq 1B,nlo)
          am = where(strmatch(result[good].r23branch,'*A*') eq 1B,nam)
          djs_oplot, mbsub[good[up]], data[good[up]].zstrong_12oh_upper, ps=8, color='red'
          djs_oplot, mbsub[good[lo]], data[good[lo]].zstrong_12oh_lower, ps=8, color='blue'
;         djs_oplot, mbsub[good[up]], data[good[up]].zstrong_12oh_lower, ps=3, color='red'
;         djs_oplot, mbsub[good[lo]], data[good[lo]].zstrong_12oh_upper, ps=3, color='blue'
          if (nam gt 1L) then djs_oplot, mbsub[good[am]], result[good[am]].zstrong_12oh, ps=8, color='green'
          sixlin, [mbsub[good[up]],mbsub[good[lo]]], [result[good[up]].zstrong_12oh,$
            result[good[lo]].zstrong_12oh], lzint, siga, lzslope, sigb
          oplot, mbaxis, poly(mbaxis,tremonti_coeff), line=2, thick=2.0
          oplot, mbaxis, poly(mbaxis,[lzint[2],lzslope[2]]), line=0, thick=2.0
          print, 'Flux/'+suffix, lzint[2], lzslope[2]
          cc = get_kbrd(1)
       endif

    endif                       ; close JUSTFLUX condition
       
; #########################
; equivalent widths
; #########################

    if keyword_set(justew) and tag_exist(data,'ZSTRONG_EW_12OH_LOWER') then begin
       
       if (not keyword_set(silent)) then splog, 'Assigning R23 branches from EWs:'

       case branchmethod of
   
          0L: begin ; just use [N II]/Ha
       
             lo = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha le n2cut) and $
               (data.zstrong_ew_12oh_lower gt -900) and (data.zstrong_ew_12oh_upper gt -900) and $
               (data.zstrong_ew_12oh_lower lt data.zstrong_ew_12oh_upper),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch_ew        = 'L_N2'
                result[lo].zstrong_ew_12oh     = data[lo].zstrong_ew_12oh_lower
                result[lo].zstrong_ew_12oh_err = data[lo].zstrong_ew_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_ew_logu     = data[lo].zstrong_ew_logu_lower
                   result[lo].zstrong_ew_logu_err = data[lo].zstrong_ew_logu_lower_err
                endif
             endif
         
             up = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha gt n2cut) and $
               (data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
               (data.zstrong_ew_12oh_upper gt data.zstrong_ew_12oh_lower),nup)
             if (nup ne 0L) then begin
                result[up].r23branch_ew        = 'U_N2'
                result[up].zstrong_ew_12oh     = data[up].zstrong_ew_12oh_upper
                result[up].zstrong_ew_12oh_err = data[up].zstrong_ew_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_ew_logu     = data[up].zstrong_ew_logu_upper
                   result[up].zstrong_ew_logu_err = data[up].zstrong_ew_logu_upper_err
                endif
             endif
             
          end
   
          1L: begin ; iterative LZ method
             
; begin by assigning branches using the O32 ratio (or put everything
; on the upper branch); on subsequent iterations, use the abundance
; that is closest to the abundance predicted by the LZ relation
; (either derived from the data, or based on Tremonti et al. 2004)

             good = where((data.zstrong_ew_12oh_upper gt -900.0) and (data.zstrong_ew_o32 gt -900.0) and (mbsub gt -900.0),ngood)
             if (ngood ne 0L) then begin
   
                mbgood = mbsub[good]
                ohgood = mbgood*0.0
                
                up = where((data[good].zstrong_ew_o32 le o32cut),nup,comp=lo,ncomp=nlo)
                if (nup ne 0L) then ohgood[up] = data[good[up]].zstrong_ew_12oh_upper
                if (nlo ne 0L) then ohgood[lo] = data[good[lo]].zstrong_ew_12oh_lower
;               ohgood = data[good].zstrong_ew_12oh_upper
                ohgood_err = ohgood*0.0

                for ii = 0L, maxiter-1L do begin
                
                   lzoh = poly(mbgood,tremonti_coeff)
;                  sixlin, mbgood, ohgood, lzint, siga, lzslope, sigb
;                  lzoh = poly(mbgood,[lzint[2],lzslope[2]])
;                  if keyword_set(debug) then begin
;                     plot, mbgood, ohgood, xsty=3, ysty=3, ps=3, $
;                       title='EW/'+strupcase(suffix), charsize=1.5, charthick=2.0, $
;                       xtitle=textoidl('M_{B}'), ytitle='12+log(O/H)'
;                     oplot, mbgood, lzoh, line=0, thick=2
;                     cc = get_kbrd(1)
;                  endif
   
                   up = where(abs(lzoh-data[good].zstrong_ew_12oh_upper) lt $
                     abs(lzoh-data[good].zstrong_ew_12oh_lower),nup,comp=lo,ncomp=nlo)
   
                   if (nup ne 0L) then begin
                      result[good[up]].r23branch_ew        = 'U_LZ'
                      result[good[up]].zstrong_ew_12oh     = data[good[up]].zstrong_ew_12oh_upper
                      result[good[up]].zstrong_ew_12oh_err = data[good[up]].zstrong_ew_12oh_upper_err
                      if keyword_set(kk04) then begin
                         result[good[up]].zstrong_ew_logu     = data[good[up]].zstrong_ew_logu_upper
                         result[good[up]].zstrong_ew_logu_err = data[good[up]].zstrong_ew_logu_upper_err
                      endif
                   endif
                   if (nlo ne 0L) then begin
                      result[good[lo]].r23branch_ew        = 'L_LZ'
                      result[good[lo]].zstrong_ew_12oh     = data[good[lo]].zstrong_ew_12oh_lower
                      result[good[lo]].zstrong_ew_12oh_err = data[good[lo]].zstrong_ew_12oh_lower_err
                      if keyword_set(kk04) then begin
                         result[good[lo]].zstrong_ew_logu     = data[good[lo]].zstrong_ew_logu_lower
                         result[good[lo]].zstrong_ew_logu_err = data[good[lo]].zstrong_ew_logu_lower_err
                      endif
                   endif
   
                   ohgood = result[good].zstrong_ew_12oh ; update the abundances
                   ohgood_err = result[good].zstrong_ew_12oh_err
   
                endfor
                   
             endif
             
          end
   
          2L: begin ; combination of [N II]/Ha and LZ assignment
   
; assign the branch according to [N II]/Ha
             
             lo = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha le n2cut) and $
               (data.zstrong_ew_12oh_lower gt -900) and (data.zstrong_ew_12oh_upper gt -900) and $
               (data.zstrong_ew_12oh_lower lt data.zstrong_ew_12oh_upper),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch_ew        = 'L_N2'
                result[lo].zstrong_ew_12oh     = data[lo].zstrong_ew_12oh_lower
                result[lo].zstrong_ew_12oh_err = data[lo].zstrong_ew_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_ew_logu     = data[lo].zstrong_ew_logu_lower
                   result[lo].zstrong_ew_logu_err = data[lo].zstrong_ew_logu_lower_err
                endif
             endif
         
             up = where((data.zstrong_niiha gt -900) and (data.zstrong_niiha gt n2cut) and $
               (data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
               (data.zstrong_ew_12oh_upper gt data.zstrong_ew_12oh_lower),nup)
             if (nup ne 0L) then begin
                result[up].r23branch_ew        = 'U_N2'
                result[up].zstrong_ew_12oh     = data[up].zstrong_ew_12oh_upper
                result[up].zstrong_ew_12oh_err = data[up].zstrong_ew_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_ew_logu     = data[up].zstrong_ew_logu_upper
                   result[up].zstrong_ew_logu_err = data[up].zstrong_ew_logu_upper_err
                endif
             endif

; fit a bisector to the data (or use the Tremonti LZ relation)          
             
             good = where((result.zstrong_ew_12oh gt -900.0) and (mbsub gt -900.0),ngood)
             if (ngood ne 0L) then begin
;               sixlin, mbsub[good], result[good].zstrong_ew_12oh, lzint, siga, lzslope, sigb
;               coeff = [lzint[2],lzslope[2]]
                coeff = tremonti_coeff
             endif 
   
; assign branches to the ambigious objects          
   
;            need = lindgen(nsubindx) & nneed = nsubindx
;            need = where((result.zstrong_ew_12oh lt -900.0) and (data.zstrong_ew_12oh_upper gt -900) and $
;              (data.zstrong_ew_12oh_lower gt -900),nneed)
             need = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900),nneed)
             if (nneed ne 0L) then begin
   
                lzoh = poly(mbsub[need],coeff)
                up = where((abs(lzoh-data[need].zstrong_ew_12oh_upper) lt abs(lzoh-data[need].zstrong_ew_12oh_lower)),$
                  nup,comp=lo,ncomp=nlo)
   
                if (nup ne 0L) then begin
                   result[need[up]].r23branch_ew        = 'U_LZ'
                   result[need[up]].zstrong_ew_12oh     = data[need[up]].zstrong_ew_12oh_upper
                   result[need[up]].zstrong_ew_12oh_err = data[need[up]].zstrong_ew_12oh_upper_err
                   if keyword_set(kk04) then begin
                      result[need[up]].zstrong_ew_logu     = data[need[up]].zstrong_ew_logu_upper
                      result[need[up]].zstrong_ew_logu_err = data[need[up]].zstrong_ew_logu_upper_err
                   endif
                endif
                if (nlo ne 0L) then begin
                   result[need[lo]].r23branch_ew        = 'L_LZ'
                   result[need[lo]].zstrong_ew_12oh     = data[need[lo]].zstrong_ew_12oh_lower
                   result[need[lo]].zstrong_ew_12oh_err = data[need[lo]].zstrong_ew_12oh_lower_err
                   if keyword_set(kk04) then begin
                      result[need[lo]].zstrong_ew_logu     = data[need[lo]].zstrong_ew_logu_lower
                      result[need[lo]].zstrong_ew_logu_err = data[need[lo]].zstrong_ew_logu_lower_err
                   endif
                endif
   
             endif
                
          end 

          3L: splog, 'BRANCHMETHOD = 3 with EWs not supported yet.'
          
          4L: begin             ; put them all on the upper branch

             good = where(data.zstrong_ew_12oh_upper gt -900.0,ngood)
             if (ngood ne 0L) then begin
                result[good].r23branch_ew        = 'U'
                result[good].zstrong_ew_12oh     = data[good].zstrong_ew_12oh_upper
                result[good].zstrong_ew_12oh_err = data[good].zstrong_ew_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[good].zstrong_ew_logu     = data[good].zstrong_ew_logu_upper
                   result[good].zstrong_ew_logu_err = data[good].zstrong_ew_logu_upper_err
                endif
             endif

          end

          5L: begin ; assign them according to R23BRANCH

             up = where((data.zstrong_ew_12oh_upper gt -900.0) and (r23branchsub eq 'U'),nup)
             if (nup ne 0L) then begin
                result[up].r23branch_ew        = 'U'
                result[up].zstrong_ew_12oh     = data[up].zstrong_ew_12oh_upper
                result[up].zstrong_ew_12oh_err = data[up].zstrong_ew_12oh_upper_err
                if keyword_set(kk04) then begin
                   result[up].zstrong_ew_logu     = data[up].zstrong_ew_logu_upper
                   result[up].zstrong_ew_logu_err = data[up].zstrong_ew_logu_upper_err
                endif
             endif

             lo = where((data.zstrong_ew_12oh_lower gt -900.0) and (r23branchsub eq 'L'),nlo)
             if (nlo ne 0L) then begin
                result[lo].r23branch_ew        = 'L'
                result[lo].zstrong_ew_12oh     = data[lo].zstrong_ew_12oh_lower
                result[lo].zstrong_ew_12oh_err = data[lo].zstrong_ew_12oh_lower_err
                if keyword_set(kk04) then begin
                   result[lo].zstrong_ew_logu     = data[lo].zstrong_ew_logu_lower
                   result[lo].zstrong_ew_logu_err = data[lo].zstrong_ew_logu_lower_err
                endif
             endif

          end

          else: begin
             splog, 'BRANCHMETHOD value not supported!'
             return, -1L
          endelse
   
       endcase 

; adopt the "average" lower/upper abundance of objects with ambiguous
; abundances, near the turn-around region; reject objects with a
; "frac" value less than fraccut (see IM_ABUNDANCE for details)

       ambig = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
         ((data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower)) or $ ; off the R23 calibration
         ((data.zstrong_ew_12oh_upper gt data.zstrong_ew_12oh_lower) and $ ; on the R23 calibration, but ambiguous
         (data.zstrong_ew_12oh_upper-data.zstrong_ew_12oh_upper_err) lt $
         (data.zstrong_ew_12oh_lower+data.zstrong_ew_12oh_lower_err)),nambig)
;      ambig = where((data.zstrong_ew_12oh_upper gt -900.0) and (data.zstrong_ew_12oh_lower gt -900.0) and $
;        (data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower),nambig)
       if (not keyword_set(silent)) then splog, '   Ambiguous: (O/H)_lower>(O/H)_upper: '+$
         string(nambig,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngalaxy),format='(F12.1)'),2)+'%).'

; reject objects whose upper- and lower-branch abundances differ by
; more than 1-sigma, where 1-sigma is given by the width of the
; overlapping histogram (12oh_avg)
       if (nambig ne 0L) then begin
          good = where((data[ambig].zstrong_ew_12oh_avg+data[ambig].zstrong_ew_12oh_avg_err gt data[ambig].zstrong_ew_12oh_lower) and $
            (data[ambig].zstrong_ew_12oh_avg-data[ambig].zstrong_ew_12oh_avg_err lt data[ambig].zstrong_ew_12oh_upper),ngood,$
            comp=reject,ncomp=nreject)
          if (not keyword_set(silent)) then splog, '     Retain: '+$
            string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%).'
          if (ngood ne 0L) then begin
             result[ambig[good]].r23branch_ew        = 'A'
             result[ambig[good]].zstrong_ew_12oh     = data[ambig[good]].zstrong_ew_12oh_avg
             result[ambig[good]].zstrong_ew_12oh_err = data[ambig[good]].zstrong_ew_12oh_avg_err
             if keyword_set(kk04) then begin
                result[ambig[good]].zstrong_ew_logu     = data[ambig[good]].zstrong_ew_logu_avg
                result[ambig[good]].zstrong_ew_logu_err = data[ambig[good]].zstrong_ew_logu_avg_err
             endif
          endif
          if (not keyword_set(silent)) then splog, '     Reject: '+$
            string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*nreject/float(nambig),format='(F12.1)'),2)+'%).'
          if (nreject ne 0L) then begin
             result[ambig[reject]].r23branch_ew        = 'Rejected'
             result[ambig[reject]].zstrong_ew_12oh     = -999.0
             result[ambig[reject]].zstrong_ew_12oh_err = -999.0
             result[ambig[reject]].zstrong_ew_logu     = -999.0
             result[ambig[reject]].zstrong_ew_logu_err = -999.0
          endif
       endif
; OLD CODE THAT USES FRACCUT
;      if (nambig ne 0L) then begin
;         good = where((data[ambig].zstrong_ew_12oh_frac gt fraccut),ngood,comp=reject,ncomp=nreject)
;         if (ngood ne 0L) then begin
;            if (not keyword_set(silent)) then splog, '     Retain (F>'+$
;              strtrim(string(fraccut,format='(F12.1)'),2)+'): '+$
;              string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
;              strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%).'
;            result[ambig[good]].r23branch_ew        = 'A'
;            result[ambig[good]].zstrong_ew_12oh     = data[ambig[good]].zstrong_ew_12oh_avg
;            result[ambig[good]].zstrong_ew_12oh_err = data[ambig[good]].zstrong_ew_12oh_avg_err
;            if keyword_set(kk04) then begin
;               result[ambig[good]].zstrong_ew_logu     = data[ambig[good]].zstrong_ew_logu_avg
;               result[ambig[good]].zstrong_ew_logu_err = data[ambig[good]].zstrong_ew_logu_avg_err
;            endif
;         endif
;         if (nreject ne 0L) then begin
;            if (not keyword_set(silent)) then splog, '     Reject (F<'+$
;              strtrim(string(fraccut,format='(F12.1)'),2)+'): '+$
;              string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
;              strtrim(string(100.0*nreject/float(nambig),format='(F12.1)'),2)+'%).'
;            result[ambig[reject]].r23branch_ew        = 'Rejected'
;            result[ambig[reject]].zstrong_ew_12oh     = -999.0
;            result[ambig[reject]].zstrong_ew_12oh_err = -999.0
;            result[ambig[reject]].zstrong_ew_logu     = -999.0
;            result[ambig[reject]].zstrong_ew_logu_err = -999.0
;         endif
;      endif
       
; average the abundances of objects in the turn-around region    
;
;      if keyword_set(test) then begin
;         ambig = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
;           ((data.zstrong_ew_12oh_upper gt data.zstrong_ew_12oh_lower) and $ ; from above
;            ((data.zstrong_ew_12oh_upper-data.zstrong_ew_12oh_upper_err) lt $
;             (data.zstrong_ew_12oh_lower+data.zstrong_ew_12oh_lower_err))) or $
;           ((data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower) and $ ; from below
;            ((data.zstrong_ew_12oh_upper+data.zstrong_ew_12oh_upper_err) gt $
;             (data.zstrong_ew_12oh_lower-data.zstrong_ew_12oh_lower_err))),nambig)
;      endif else begin
;         ambig = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
;           ((data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower) and $ ; from below
;            ((data.zstrong_ew_12oh_upper+data.zstrong_ew_12oh_upper_err) gt $
;             (data.zstrong_ew_12oh_lower-data.zstrong_ew_12oh_lower_err))),nambig)
;      endelse
;
;      if (nambig ne 0L) then begin
;         result[ambig].r23branch_ew = 'A'
;         for i = 0L, nambig-1L do begin
;            oh     = [data[ambig[i]].zstrong_ew_12oh_upper,data[ambig[i]].zstrong_ew_12oh_lower]
;            oh_err = [data[ambig[i]].zstrong_ew_12oh_upper_err,data[ambig[i]].zstrong_ew_12oh_lower_err]
;            result[ambig[i]].zstrong_ew_12oh     = total(oh/oh_err^2)/total(1.0/oh_err^2)
;            result[ambig[i]].zstrong_ew_12oh_err = 1.0/sqrt(total(1.0/oh_err^2))
;            print, oh[0], oh[1], result[ambig[i]].zstrong_ew_12oh, result[ambig[i]].zstrong_ew_12oh_err
;         endfor
;      endif
;
;; reject objects having O/H_upper<O/H_lower that are not within
;; one-sigma of one another (these usually have R23>10) (this step
;; should not be necessary)
;    
;       reject = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
;         ((data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower) and $ 
;          ((data.zstrong_ew_12oh_upper+data.zstrong_ew_12oh_upper_err) lt $ ; note!
;           (data.zstrong_ew_12oh_lower-data.zstrong_ew_12oh_lower_err))),nreject)
;       if (nreject ne 0L) then begin
;;         splog, 'Rejecting (EW)  : (O/H)_lower>(O/H)_upper: '+string(nreject,format='(I0)')+'/'+$
;;           string(ngalaxy,format='(I0)')+' ('+string(round(100.0*nreject/float(ngalaxy)),format='(I0)')+'%).'
;          result[reject].r23branch_ew        = 'Rejected'
;          result[reject].zstrong_ew_12oh     = -999.0
;          result[reject].zstrong_ew_12oh_err = -999.0
;       endif

       if keyword_set(debug) and (n_elements(mb) ne 0L) then begin
          good = where((mbsub gt -900.0) and (result.zstrong_ew_12oh gt -900.0))
          plot, mbsub[good], result[good].zstrong_ew_12oh, ps=8, xsty=3, ysty=3, $
            title='EW/'+strupcase(suffix)+' - Final', charsize=1.5, charthick=2.0, $
            xtitle=textoidl('M_{B}'), ytitle='12+log(O/H)'
          up = where(strmatch(result[good].r23branch_ew,'*U*') eq 1B,nup) 
          lo = where(strmatch(result[good].r23branch_ew,'*L*') eq 1B,nlo)
          am = where(strmatch(result[good].r23branch_ew,'*A*') eq 1B,nam)
          djs_oplot, mbsub[good[up]], data[good[up]].zstrong_ew_12oh_upper, ps=8, color='red'
          djs_oplot, mbsub[good[lo]], data[good[lo]].zstrong_ew_12oh_lower, ps=8, color='blue'
;         djs_oplot, mbsub[good[up]], data[good[up]].zstrong_ew_12oh_lower, ps=3, color='blue'
;         djs_oplot, mbsub[good[lo]], data[good[lo]].zstrong_ew_12oh_upper, ps=3, color='red'
          if (nam gt 1L) then djs_oplot, mbsub[good[am]], result[good[am]].zstrong_ew_12oh, ps=8, color='green'
          sixlin, [mbsub[good[up]],mbsub[good[lo]]], [result[good[up]].zstrong_ew_12oh,$
            result[good[lo]].zstrong_ew_12oh], lzint, siga, lzslope, sigb
          oplot, mbaxis, poly(mbaxis,[lzint[2],lzslope[2]]), line=0, thick=2.0
          oplot, mbaxis, poly(mbaxis,tremonti_coeff), line=2, thick=2.0
          print, 'EW/'+suffix, lzint[2], lzslope[2]
          cc = get_kbrd(1)
       endif

    endif                       ; close JUSTEW condition
       
; rename the tags and return    

    allresult[subindx] = result
    allresult = im_struct_trimtags(allresult,select=tag_names(allresult),newtags=finaltags)

return, allresult
end
    

;   diff = abs(abund[index].zstrong_ew_12oh_m91_upper-abund[index].zstrong_ew_12oh_m91_lower)
;   ambig1 = where((abund[index].zstrong_ew_12oh_m91_upper lt abund[index].zstrong_ew_12oh_m91_lower) and $
;     diff lt sqrt(abund[index].zstrong_ew_12oh_m91_upper_err^2.0+abund[index].zstrong_ew_12oh_m91_lower_err^2.0))
;   ambig2 = where((abund[index].zstrong_ew_12oh_m91_upper lt abund[index].zstrong_ew_12oh_m91_lower) and $
;     diff gt sqrt(abund[index].zstrong_ew_12oh_m91_upper_err^2.0+abund[index].zstrong_ew_12oh_m91_lower_err^2.0))

;   ambig = where((abund[index].zstrong_ew_12oh_m91_upper lt abund[index].zstrong_ew_12oh_m91_lower) and $
;     (abund[index].zstrong_ew_12oh_m91_frac gt 0.2))
;   plot, abund[index[ambig]].zstrong_ew_r23, abund[index[ambig]].zstrong_ew_12oh_m91_avg, $
;     ps=4, xsty=3, ysty=3, xr=[0.6,1.8], yr=[7.8,9.2]
;   ploterror, abund[index[ambig]].zstrong_ew_r23, abund[index[ambig]].zstrong_ew_12oh_m91_avg, $
;     abund[index[ambig]].zstrong_ew_r23_err, abund[index[ambig]].zstrong_ew_12oh_m91_avg_err, $
;     ps=4, xsty=3, ysty=3, xr=[0.6,1.8], yr=[7.8,9.2]
;   ploterror, abund[index[ambig]].zstrong_ew_r23, abund[index[ambig]].zstrong_ew_12oh_m91_avg, $
;     abund[index[ambig]].zstrong_ew_12oh_m91_avg_err, ps=4, xsty=3, ysty=3, xr=[0.6,1.8], yr=[7.8,9.2]
    
;;; ---------------------------------------------------------------------------    
;;    
;;    good = where((data.zstrong_12oh_upper gt -900.0) and (data.zstrong_12oh_lower gt -900.0) and (mbsub gt -900.0),ngood)
;;    mbgood = mbsub[good]
;;    djs_plot, data[good].zstrong_r23, data[good].zstrong_12oh_upper, xsty=3, ysty=3, ps=3, color='red', yrange=[7,9.5]
;;    djs_oplot, data[good].zstrong_r23, data[good].zstrong_12oh_lower, ps=3, color='blue'
;;
;;; average the abundances of objects in the turn-around region    
;;
;;    if keyword_set(test) then begin
;;       ambig = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
;;         ((data.zstrong_12oh_upper gt data.zstrong_12oh_lower) and $ ; from above
;;          ((data.zstrong_12oh_upper-data.zstrong_12oh_upper_err) lt $
;;           (data.zstrong_12oh_lower+data.zstrong_12oh_lower_err))) or $
;;         ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower) and $ ; from below
;;          ((data.zstrong_12oh_upper+data.zstrong_12oh_upper_err) gt $
;;           (data.zstrong_12oh_lower-data.zstrong_12oh_lower_err))),nambig)
;;    endif else begin
;;       ambig = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
;;         ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower) and $ ; from below
;;          ((data.zstrong_12oh_upper+data.zstrong_12oh_upper_err) gt $
;;           (data.zstrong_12oh_lower-data.zstrong_12oh_lower_err))),nambig)
;;    endelse
;;    if (nambig ne 0L) then begin
;;       result[ambig].r23branch = 'A'
;;       for i = 0L, nambig-1L do begin
;;          oh     = [data[ambig[i]].zstrong_12oh_upper,data[ambig[i]].zstrong_12oh_lower]
;;          oh_err = [data[ambig[i]].zstrong_12oh_upper_err,data[ambig[i]].zstrong_12oh_lower_err]
;;          result[ambig[i]].zstrong_12oh     = total(oh/oh_err^2)/total(1.0/oh_err^2)
;;          result[ambig[i]].zstrong_12oh_err = 1.0/sqrt(total(1.0/oh_err^2))
;;       endfor
;;       djs_oplot, data[ambig].zstrong_r23, result[ambig].zstrong_12oh, ps=4, syms=0.5, color='yellow'
;;       djs_oplot, data[ambig].zstrong_r23, data[ambig].zstrong_12oh_upper, ps=3, color='yellow'
;;       djs_oplot, data[ambig].zstrong_r23, data[ambig].zstrong_12oh_lower, ps=3, color='yellow'
;;    endif
;;
;;    reject = where((result.zstrong_12oh gt -900.0) and $
;;      (data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900.0) and $
;;      ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower) and $ 
;;       ((data.zstrong_12oh_upper+data.zstrong_12oh_upper_err) lt $ ; note!
;;        (data.zstrong_12oh_lower-data.zstrong_12oh_lower_err))),nreject)
;;    if (nreject ne 0L) then begin
;;       djs_oplot, data[reject].zstrong_r23, result[reject].zstrong_12oh, ps=4, color='purple'
;;       result[reject].r23branch        = 'Rejected'
;;       result[reject].zstrong_12oh     = -999.0
;;       result[reject].zstrong_12oh_err = -999.0
;;    endif
;;
;;; now compare the abundances to the SDSS LZ relation
;;
;;    djs_plot, mbgood, data[good].zstrong_12oh_upper, xsty=3, ysty=3, ps=3, color='red', yrange=[7,9.5]
;;    djs_oplot, mbgood, data[good].zstrong_12oh_lower, ps=3, color='blue'
;;    lzoh = poly(mbgood,tremonti_coeff)
;;    up = where((abs(lzoh-data[good].zstrong_12oh_upper) lt abs(lzoh-data[good].zstrong_12oh_lower)),nup,comp=lo,ncomp=nlo)
;;    djs_oplot, mbgood[up], data[good[up]].zstrong_12oh_upper, ps=3, color='yellow'
;;    djs_oplot, mbgood[lo], data[good[lo]].zstrong_12oh_lower, ps=3, color='green'
;;    
;; ---------------------------------------------------------------------------
    
