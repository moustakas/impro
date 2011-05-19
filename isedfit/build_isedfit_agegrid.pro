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
  minage=minage, maxage=maxage, lookback=lookback, linear=linear, $
  debug=debug

    if keyword_set(linear) then log = 0 else log = 1
    
    if (n_elements(inage) eq 0L) then begin
       if (n_elements(minage) eq 0) then minage = 0.01D
       if (n_elements(maxage) eq 0) then maxage = 13.5D
       if (n_elements(nage) eq 0) then nage = 100
       inage = range(minage,maxage,nage,log=log)
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
    if (nsamp gt 5) then begin                ; if NAGE is really small....
       fact = range(0.1,4.0,(nsamp-1)/2,/log) ; 0.1-4-sigma

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
