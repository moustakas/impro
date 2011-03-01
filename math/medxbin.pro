;+
; NAME:
;       MEDXBIN
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;
;-

pro medxbin, x, y, binsz, medx = medx, medy = medy, sigy = sigy, $
  clip = clip, niter = niter, minmax = minmax, thresh = thresh, $
  distr = distr, minpoints = minpoints

    nbins = round((max(x) - min(x)) / binsz) + 1L
    medx = fltarr(nbins) + 'NaN'
    medy = fltarr(nbins) + 'NaN'
    sigy = fltarr(nbins) + 'NaN'
    distr = fltarr(4, nbins) + 'NaN'

    if (n_elements(thresh) eq 0L) then thresh = 15
    if (n_elements(minpoints) eq 0L) then minpoints = 10L
    if not keyword_set(niter) then niter = 5
    if not keyword_set(clip) then clip = 5

    for i = 0, nbins - 1 do begin

       binctr = i * binsz + min(x) ;+ (binsz/2.0)

       inbin = where((x gt binctr - (binsz/2.0)) and (x le binctr + (binsz/2.0)),ninbin)

       if (ninbin ne 0L) then begin

          if n_elements(inbin) lt minpoints then continue

          if keyword_set(minmax) then begin
             clip = where(y[inbin] ge minmax[0] and y[inbin] le minmax[1])
             if clip[0] ne -1 then inbin = inbin[clip]
          endif

          ok = where(finite(y[inbin]), nok)
;         print, nok, thresh

          if nok lt thresh then continue

          inbin = inbin[ok]
          medx[i] = median(x[inbin])
          medy[i] = median(y[inbin])

;         meanclip, y[inbin], mean, sigyi, clip=clip, maxiter=niter, converge=0.001, verbose=0
;         sigyi = robust_sigma(y[inbin])
;         sigy[i] = sigyi

          yhist = y[inbin]
          yhist = yhist[sort(yhist)]
          linterp, findgen(nok)/nok, yhist, [0.025, 0.16, 0.84, 0.975], distri
          distr[*,i] = distri
          sigy[i] = (distri[2] - distri[1]) / 2

       endif
          
       print, i, medx[i], medy[i]

    endfor

    ok = where(finite(medy) eq 1)
    if ok[0] eq -1 then begin
       medx = [0, 0]
       medy = [0, 0]
       sigy = [0, 0]
    endif else begin
       medx = medx[ok]
       medy = medy[ok]
       sigy = sigy[ok]
       distr = distr[*,ok]
    endelse

return    
end

