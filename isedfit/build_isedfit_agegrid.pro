;+
; NAME:
;   BUILD_ISEDFIT_AGEGRID()
;
; PURPOSE:
;   Build an age (time) grid that deals with bursts correctly. 
;
; INPUTS: 
;   info - iSEDfit-style star formation history structure 
;     tau - characteristic e-folding time [Gyr]
;     nburst - number of bursts
;     tburst - time each burst begins [Gyr]
;     dtburst - burst duration [Gyr]
;     fburst - mass fraction of each burst
;     truncate - construct an age vector for a truncated burst (only
;       used for NBURST=1)
;
; OPTIONAL INPUTS: 
;   inage - input age vector to use as a starting point [Gyr]; note
;     that you either have to pass INAGE *or* NAGE, MINAGE, MAXAGE
;     (with INAGE dominating) 
;   nage - number of elements in OUTAGE (default 100)
;   minage - minimum age [Gyr] (default 0.05)
;   maxage - maximum age [Gyr] (default 13.0)
;
; KEYWORD PARAMETERS: 
;   lookback - convert OUTAGE to a lookback time
;   debug - make a simple debugging plot and wait for a keystroke 
;
; OUTPUTS: 
;   outage - output age vector [Gyr]
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Jan 01, UCSD
;
; Copyright (C) 2011, John Moustakas
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

function build_isedfit_agegrid, info, inage=inage, nage=nage, $
  minage=minage, maxage=maxage, lookback=lookback, debug=debug

    if (n_elements(inage) eq 0L) then begin
       if (n_elements(minage) eq 0) then minage = 0.05D
       if (n_elements(maxage) eq 0) then maxage = 13.0D
       if (n_elements(nage) eq 0) then nage = 100
       inage = range(minage,maxage,nage,/log)
       if keyword_set(lookback) then inage = $
         reverse(maxage-(inage-min(inage)))
    endif
    minage = min(inage)
    maxage = max(inage)
    nage = n_elements(inage)
    
; some convenient variables    
    nb = info.nburst
    if (nb eq 0) then return, inage

    tb = info.tburst[0:nb-1]
    dtb = info.dtburst[0:nb-1]

; demand that at least NSAMP points are used to sample the burst, but
; (in the limit of lots of bursts) no more than 80% of NAGE
    nsamp = long(round(nage*0.12)<(0.8*nage/float(nb))) ; ages per burst
    nsamp = nsamp+1*(odd(nsamp) eq 0) ; ensure odd
    if (nsamp gt 5) then begin ; if NAGE is really small....
       fact = range(0.3,4.0,(nsamp-1)/2,/log) ; 0.3-4-sigma

       for ib = 0, nb-1 do begin
          burstage1 = [reverse(tb[ib]-dtb[ib]*fact),tb[ib],tb[ib]+dtb[ib]*fact]
          exclude = where((burstage1 le minage) or (burstage1 ge maxage),$
            nexclude,comp=good,ncomp=ngood)
          if (ngood ne 0) then begin
             burstage1 = burstage1[good] ; NGOOD=0 means the burst has not happened yet
             if (burstage1[0] ne -1) then begin
                if (n_elements(burstage) eq 0) then burstage = burstage1 else $
                  burstage = [burstage,burstage1]
             endif
          endif
       endfor

; now get rid of NBURSTAGE random ages from the old age vector
       ntoss = n_elements(burstage)
       if (ntoss ne 0) then begin
          these = cmset_op(inage,'and',/not2,burstage,/index) ; remove TB from INAGE
          nthese = n_elements(these)
          keep = random_indices(nthese,nage-ntoss)
          outage = [inage[these[keep]],burstage]
       endif else outage = inage
       outage = outage[sort(outage)]
    endif else outage = inage

    if (n_elements(outage) ne nage) then message, 'Bad bad bad!'
    
; truncate the last burst?  require at least 5 age samplings of the
; exponential tail 
    if tag_exist(info,'tauburst') then begin
       if (info.tauburst gt 0.0) then begin
          post = where(outage ge tb[nb-1],npost,comp=pre,ncomp=npre)
          if (npost gt 0) then begin
             npad = 5
             padage = tb[nb-1]-info.tauburst*alog((1D0-(range(0.1D,1D,npad,/log)-0.1D)))
             keep = where(padage lt maxage,nkeep) ; could break if TB=MAXAGE
             if (nkeep gt 0) then begin
                padage = padage[keep]
                these = where(outage gt max(padage),nthese)
                if (nthese gt 0) then padage = [padage,outage[these]]
;  get rid of the "spare" ages from before the exponential
                if (n_elements(padage) gt npost) then begin
                   ran = random_indices(npre,npre-(n_elements(padage)-npost))
                   outage = [outage[pre[ran]],padage]
                   outage = outage[sort(outage)]
                endif
             endif
          endif
       endif 
    endif 
    
    if keyword_set(debug) then begin
       plot, outage, inage, ps=-6, /xlog, /ylog, yr=minmax(inage), xsty=3, ysty=3
       for kk = 0, nb-1 do djs_oplot, tb[kk]*[1,1], 10^!y.crange, color='red'
       cc = get_kbrd(1)
    endif

    uu = uniq(outage,sort(outage))
    if (n_elements(uu) lt n_elements(outage)) then message, 'Bad bad'
    
return, outage
end

;   if ((max(tb+dtb)-maxage) gt 1D-4) then begin
;      maxage = max(tb+dtb)+pad
;      splog, 'Warning: TB+DTB exceeds MAXAGE! ', maxage
;      inage = range(minage,maxage,nage,/log)
;      if keyword_set(lookback) then inage = $
;        reverse(maxage-(inage-min(inage)))
;   endif
;
;; compute the number of pixels to assign between bursts based on the
;; fractional time; require at least *five* pixels between bursts
;    tb1 = fltarr(2*nb+1)-1
;    tb2 = tb1
;    tb1[0] = minage
;    tb2[0] = tb[0]-pad
;    cc = 0
;    for ii = 1, 2*nb do begin
;       if odd(ii) then begin
;          tb1[ii] = tb[cc]
;          tb2[ii] = tb[cc]+dtb[cc]
;       endif else begin
;; check for overlapping bursts
;          if ((tb[cc]+dtb[cc]) lt tb[cc+1]) then $
;            tb1[ii] = tb[cc]+dtb[cc]+pad
;          if (ii eq 2*nb) then tb2[ii] = maxage else tb2[ii] = tb[cc+1]-pad
;          cc++
;       endelse
;    endfor
;
;    tb1 = [minage,tb,tb+dtb+pad]
;    tb2 = [tb-pad,tb+dtb,maxage]
;;   tb1 = tb1[sort(tb1)]
;;   tb2 = tb2[sort(tb2)]
;    npix = round(nage*(tb2-tb1)/(maxage-minage))>3
;    
;; deal with "extra" (or too few) pixel(s)
;    diff = long(nage-total(npix,/double))
;    if (abs(diff) gt 0) then junk = max(npix,this) else $
;      junk = min(npix,this)
;    npix[this] = npix[this] + diff
;    pix = [0L,long(total(npix,/cumu))] ; bin boundaries
;;   niceprint, tb1, tb2, npix
;
;    neg = where(npix le 0)
;    if (neg[0] ne -1) then message, 'This is bad'
;    
;; build the output time grid; note that we ensure integer pixels
;; corresponding to the beginning and end of each burst
;    outage = inage*0D
;    for jj = 0, 2*nb do begin
;       thisage = range(tb1[jj],tb2[jj],npix[jj],/log) ; log-spacing
;       if keyword_set(lookback) then outage[pix[jj]:pix[jj+1]-1] = $
;         reverse(max(thisage)-(thisage-min(thisage))) else $
;           outage[pix[jj]:pix[jj+1]-1] = thisage    
;    endfor
;    outage = outage[sort(outage)] ; resort the output
;;   niceprint, findex(outage,tb), findex(outage,tb+dtb)


;;
;;
;;       burstage = dblarr(nsamp*nb)
;;       for ii = 0, nb-1 do burstage[ii*nsamp:nsamp*(ii+1)-1] = $
;;         [reverse(tb[ii]-dtb[ii]*fact),tb[ii],tb[ii]+dtb[ii]*fact]
;;       exclude = where((burstage lt minage) or (burstage ge maxage),$
;;         nexclude,comp=good,ncomp=ngood)
;;
;;       if (ngood eq 0) then return, inage ; NGOOD=0 means the burst has not happened yet
;;       burstage = burstage[good]
;;
;;; throw out ages from the old (input) age vector, unless we throw out
;;; more than we're adding (usually true for a very early burst) 
;;       keep = where((inage lt min(burstage)) or (inage gt max(burstage)),nkeep)
;;       if (nkeep gt 0) and (nkeep lt ngood) then outage = [inage[keep],burstage] else $
;;         outage = [inage,burstage]
;;       outage = outage[uniq(outage,sort(outage))]
;;;      outage = outage[sort(outage)]
;;       isburst = fix(outage*0)
;;       for bb = 0, n_elements(outage)-1 do isburst[bb] = total(outage[bb] eq burstage) ge 1
;;    
;;; at this point our age vector contains (NSAMP*NB)<NGOOD too many
;;; ages, so remove the "extra" ages from periods when the "galaxy"
;;; isn't bursting, but do not touch MINAGE and MAXAGE; override
;;; truncated bursts, if any 
;;       junk = isedfit_reconstruct_sfh(info,age=outage,$
;;         sfhburst=sfhburst,/notruncate)
;;       sfhburst = sfhburst/max(sfhburst)
;;       donottouch = (outage eq minage) or (outage eq maxage) or (isburst eq 1)
;;       calm = where((sfhburst lt 0.01) and (donottouch eq 0),ncalm,$
;;         comp=bursty,ncomp=nbursty)
;;       if (ngood gt ncalm) then message, 'Your star formation history is too bursty!'
;;
;;      djs_plot, outage, sfhburst, psym=-6, /xlog, ysty=3, xsty=3, xr=tb[0]+[-2,2]
;;      djs_oplot, outage[calm], sfhburst[calm], psym=6, color='cyan'
;;
;;       keep = random_indices(ncalm,ncalm-ngood)
;;       outage = [outage[calm[keep]],outage[bursty]]
;;       outage = outage[sort(outage)]
;;    endif else outage = inage
;;    if (n_elements(outage) ne nage) then message, 'Bad bad bad!'
;;    
