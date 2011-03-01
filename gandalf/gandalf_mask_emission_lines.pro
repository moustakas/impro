function gandalf_mask_emission_lines, npix, Vsys, linepars, $
  velscale, l0_gal, lstep_gal, sigma=sigma, l_rf_range=l_rf_range
; Return a list of goodpixels to fit that excludes regions potentially
; affected by gas emission and by sky lines. Unless the log10 keyword
; is specified, wavelength values are assumed to be ln-rebinned, and
; are defined by the l0_gal, lstep_gal, npix parameters. The position of
; gas and sky emission lines is set by the input linepars
; structure and the width of the mask by the sigma parameter. If a
; sigma value is not passed than the width of each line is taken from
; the linepars structure.
;
; The rest-frame fitting wavelength range can be manually restricted
; using the l_rf_range keyword to pass min and max observed
; wavelength. Typically used to exclude regions at either side of
; spectra.

    nline = n_elements(linepars.i)
    
; speed of light
    c = 299792.458d
; define good pixels array
    goodpixels = range(0,npix-1) 
; if set, exclude regions at either ends of the spectra using the keyword l_rf_range
    if keyword_set(l_rf_range) then begin
       pix0 = ceil((alog(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
       pix1 = ceil((alog(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
       goodpixels = range(max([pix0,0]),min([pix1,npix-1])) ; NEW - V1.3 
    endif

    tmppixels  = goodpixels

; looping over the listed emission-lines and mask those tagged with an
; 'm' for mask. Mask sky lines at rest-frame wavelength
    for i = 0, nline-1 do begin
       if (linepars[i].action eq 'm') then begin
;         print,'--> masking ' + linepars.name[i]
          if (linepars[i].name ne 'sky') then $
            meml_cpix = ceil((alog(linepars[i].lambda)-l0_gal)/lstep_gal+Vsys/velscale)
; sky lines are at rest-frame
          if (linepars[i].name eq 'sky') then $
            meml_cpix = ceil((alog(linepars[i].lambda)-l0_gal)/lstep_gal) 
; set the width of the mask in pixels using either 3 times the sigma
; of each line in the emission-line setup or the provided sigma value;
; mask broad lines by a larger amount
          if strmatch(linepars[i].name,'*broad*',/fold) then $
            factor = 20.0 else factor = 3.0 ; 10-sigma vs 3-sigma
;         factor = 3.0
          if keyword_set(sigma) then $
            msigma = factor*sigma/velscale else $
            msigma = factor*linepars[i].s/velscale
          meml_bpix = meml_cpix - msigma
          meml_rpix = meml_cpix + msigma
          w = where(goodpixels ge meml_bpix and goodpixels le meml_rpix) 
          if (w[0] ne -1) then begin
             tmppixels[w] = -1 
          endif else begin 
;            print, linepars[i].name+' is outside your wavelength range. We shall ignore it' ; NEW - V1.3
             linepars[i].action = 'i'                                         ; NEW - V1.3 
          endelse
       endif
    endfor
    w = where(tmppixels ne -1)
    goodpixels = goodpixels[w]

; check for lines that are being fitted that are outside the
; wavelength range, and ignore those, too
    for i = 0, nline-1 do begin
       if (linepars[i].action eq 'f') then begin
          if (linepars[i].name ne 'sky') then $
            meml_cpix = ceil((alog(linepars[i].lambda)-l0_gal)/lstep_gal+Vsys/velscale)
; set the width of the mask in pixels using either 3 times the sigma
; of each line in the emission-line setup or the provided sigma value;
; mask broad lines by a larger amount
          if strmatch(linepars[i].name,'*broad*',/fold) then $
            factor = 20.0 else factor = 3.0 ; 10-sigma vs 3-sigma
;         factor = 3.0
          if keyword_set(sigma) then $
            msigma = factor*sigma/velscale else $
            msigma = factor*linepars[i].s/velscale
          meml_bpix = meml_cpix - msigma
          meml_rpix = meml_cpix + msigma
          if (meml_bpix le 2.0) or (meml_rpix ge (npix-2.0)) then begin
;            print, linepars[i].name+' is outside your wavelength range. We shall ignore it' ; NEW - V1.3
             linepars[i].action = 'i' ; NEW - V1.3 
          endif             
       endif 
    endfor

; if all the lines have been dropped then we're done
    if (total(linepars.action eq 'i') eq nline) then $
      return, goodpixels

; now we have to loop back through and check to see if any of the blue
; rest-wavelength Balmer or forbidden lines are tied to a line that
; has just been dropped (e.g., H-gamma and H-alpha tied to H-beta, or
; [NII] 6548 tied to [NII] 6584)
    sat = where(strmatch(linepars.kind,'*d*'),nsat)
    for ii = 0, nsat-1 do begin
       big = where(linepars.i eq strtrim(strmid(linepars[sat[ii]].kind,1),2),nbig)
       if (strtrim(linepars[big].action,2) eq 'i') then linepars[sat[ii]].action = 'i'
    endfor

; narrow lines
    tied = where((strmatch(linepars.fit,'*t*') eq 1) and $
      (strmatch(linepars.name,'*broad*',/fold) eq 0),ntied)
    for ii = 0, ntied-1 do begin
       this = where(linepars.i eq strtrim(strmid(linepars[tied[ii]].fit,1),2),nthis)
       if (nthis ne 0) then begin
; if [OIII] is tied then choose the bluest rest-frame line that is
; *not* tied
          if (linepars[this].action eq 'i') then begin
             these = where($
               ((strtrim(linepars.action,2) eq 'f') or $
               (strtrim(linepars.action,2) eq 'm')) and $
               (strmatch(linepars.kind,'d*') eq 0),nthese)
             if (nthese eq 0L) then message, 'Fix me!'
             bluest = min(linepars[these].lambda,blueindx)
             linepars[tied[ii]].fit = 't'+strtrim(linepars[these[blueindx]].i,2)
             linepars[these[blueindx]].fit = 'f'
;; first try tying to [OII]
;             oii = where(strmatch(linepars.name,'*oii_*',/fold),noii)
;             if (noii ne 0) then begin
;                linepars[tied[ii]].fit = 't'+strtrim(linepars[oii].i,2)
;                linepars[oii].fit = 'f'
;             endif
          endif
       endif
    endfor

; broad lines; note that the broad lines remain tied even on the
; second iteration, except for Mg II
    tied = where((strmatch(linepars.fit,'*t*') eq 1) and $
      (strmatch(linepars.name,'*broad*',/fold) eq 1),ntied)
    for ii = 0, ntied-1 do begin
       this = where(linepars.i eq strtrim(strmid(linepars[tied[ii]].fit,1),2),nthis)
       if (nthis ne 0) then begin
          if (linepars[this].action eq 'i') then begin
; very brittle code!!
; Balmer lines
             hb = where(strmatch(linepars.name,'*h_beta_broad*',/fold),nhb)
             hg = where(strmatch(linepars.name,'*h_gamma_broad*',/fold),nhg)
             if (linepars[hb].action eq 'i') then begin
                linepars[tied[ii]].fit = 't'+strtrim(linepars[hg].i,2)
                linepars[hg].fit = 'f'
                linepars[tied[ii]].fit_iter2 = 't'+strtrim(linepars[hg].i,2)
                linepars[hg].fit_iter2 = 'f'
             endif else begin
                linepars[tied[ii]].fit = 't'+strtrim(linepars[hb].i,2)
                linepars[hb].fit = 'f'
                linepars[tied[ii]].fit_iter2 = 't'+strtrim(linepars[hb].i,2)
                linepars[hb].fit_iter2 = 'f'
             endelse
; Mg II
             if strmatch(linepars[tied[ii]].name,'*mgii*',/fold) then begin
                linepars[tied[ii]].fit_iter2 = 'f'
                if (linepars[hb].action eq 'i') then begin
;                  linepars[tied[ii]].fit = 'f'
                   linepars[tied[ii]].fit = 't'+strtrim(linepars[hg].i,2)
                endif
             endif
          endif
       endif
    endfor

;   for ii = 0, n_elements(linepars.name)-1 do begin
;      istied = where((strtrim(linepars[ii].action,2) eq 'f') and $
;        ((strmid(strtrim(linepars[ii].kind,2),0,1) eq 'd') or $
;        (strmid(strtrim(linepars[ii].fit,2),0,1) eq 't')),nistied)
;      if (nistied ne 0) then begin
;         this = where((linepars[istied].action eq 'f'),nthis)
;         if (nthis ne 0) then begin
;            this = this[nthis-1] ; pick the reddest not-ignored line
;            linepars[istied].fit = 't'+strtrim(linepars[istied[this]].i,2)
;         endif
;      endif
;   endfor

return, goodpixels
end

