;+
; NAME:
;   IM_PLOTHIST
;
; PURPOSE:
;   Generate a histogram plot with optional weights.
;
; INPUTS: 
;   arr - input array [NPTS]
;
; OPTIONAL INPUTS: 
;   weight - weight for each input point [NPTS]
;   bin - histogram binsize
;   edge - desired centering: 0=center, -1=left, +1=right (default 0)
;   psym - plot symbol (default 10)
;   normfactor - normalize the histogram by this factor
;   fcolor - fill color (see PLOTHIST)
;   fline - fill linestyle (see PLOTHIST)
;   fspacing - fill line spacing (see PLOTHIST)
;   fpattern - (see PLOTHIST)
;   forientation - line orientation anlge(see PLOTHIST)
;   extra - inputs for IM_HIST1D)
;
; KEYWORD PARAMETERS: 
;    cumulative - generate the cumulative distribution function
;    overplot - render on an existing plot
;    fill - use polyfill to fill in the histogram
;    noplot - just build the histograms; do not plot
;    fraction - normalize the histogram by the total area
;
; OUTPUTS: 
;   xhist - binned x-axis
;   yhist - binned y-axis
;
; COMMENTS:
;   Combines the best of HOGG_PLOTHIST and PLOTHIST (in my opinion!) 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;    J. Moustakas, 2007 Apr 02, NYU
;
; Copyright (C) 2007, John Moustakas
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

pro im_plothist, arr, xhist, yhist, bin=bin1, edge=edge, weight=weight, $
  psym=psym, cumulative=cumulative, normfactor=normfactor, noplot=noplot, $
  overplot=overplot, fill=fill, fcolor=fcolor, fline=fline, fspacing=fspacing, $
  fpattern=fpattern, forientation=forientation, fraction=fraction, peak=peak, $
  locations=locations, xlog=xlog, logbins=logbins, _extra=extra

    ndata = n_elements(arr)
    if (ndata eq 0L) then begin
       doc_library, 'im_plothist'
       return
    endif
    
    dtype = size(arr,/type)
    if (n_elements(weight) eq 0L) then weight = dblarr(ndata)+1.0

    if (n_elements(bin1) eq 0) then begin
       if keyword_set(logbins) then $
         bin = (max(alog10(arr))-min(alog10(arr)))/double(ceil(0.3D*sqrt(ndata))) else $
           bin = (max(arr)-min(arr))/double(ceil(0.3D*sqrt(ndata)))
       if (dtype eq 4) then bin = float(bin)
    endif else begin
       bin = float(abs(bin1))
    endelse

; if ARR is an integer and BIN is floating point then histogram()
; barfs; deal with that here
    atype = size(arr,/type)
    btype = size(bin,/type)
    if ((atype eq 2) or (atype eq 3)) and ((btype eq 4) or (btype eq 5)) then $
      amult = 1.0 else amult = 1

;   if (size(bin,/type) ne size(arr,/type)) then begin
;      message, 'BIN and ARR datatypes must match!'
;   endif

; allow for logarithmic bins (code contributed by J. Aird)
    if keyword_set(logbins) then begin
       xlog = 1
       yhist1 = im_hist1d(alog10(amult*arr),weight,binsize=bin,obin=xhist1,$
         binedge=edge,_extra=extra,locations=loglocations)
       locations = 10.0^loglocations

       xhist = reform(transpose(rebin(locations,n_elements(locations),2)),n_elements(locations)*2)
       yhist = shift(reform(transpose(rebin(yhist1,n_elements(yhist1),2)),n_elements(yhist1)*2),+1)
       yhist[0] = 0
       yhist[n_elements(yhist)-1] = 0
    endif else begin
       yhist = im_hist1d(amult*arr,weight,binsize=bin,obin=xhist,$
         binedge=edge,_extra=extra)
    endelse
    n_hist = n_elements(yhist)
    
    if keyword_set(fraction) then normfactor = total(yhist)
    if (n_elements(normfactor) ne 0L) then yhist = yhist/float(normfactor)
    if keyword_set(peak) then yhist = yhist/max(yhist)

    if keyword_set(noplot) then return

    if (n_elements(extra) ne 0) then begin
       if tag_exist(extra,'xtitle') then extra.xtitle = textoidl(extra.xtitle)
       if tag_exist(extra,'ytitle') then extra.ytitle = textoidl(extra.ytitle)
       if tag_exist(extra,'color') then extra.color = im_color(extra.color)
    endif
    if (n_elements(psym) eq 0L) then psym = 10

    if keyword_set(overplot) then begin
       if keyword_set(cumulative) then begin
          yhist = total(yhist,/cumulative)/total(yhist)
          oplot, xhist, yhist, _extra=extra 
       endif else begin
          oplot, [xhist[0]-bin,xhist,xhist[n_hist-1]+bin], [0,yhist,0], $ 
            psym=psym, _extra=extra
       endelse
    endif else begin
       if keyword_set(cumulative) then begin
          yhist = total(yhist,/cumulative)/total(yhist)
          plot, xhist, yhist, xlog=xlog, _extra=extra 
       endif else begin
          plot, [xhist[0]-bin,xhist,xhist[n_hist-1]+bin], [0,yhist,0],  $ 
            psym=psym, xlog=xlog, _extra=extra
       endelse
    endelse

    if keyword_set(fill) then begin
       if keyword_set(logbins) then begin
;      if keyword_set(logbins) and n_elements(bin1) eq 0 then begin
          xfill = [xhist[0],xhist,xhist[n_hist-1]]
          yfill = [0.0,yhist,0.0]
;         xfill = transpose([[Xhist-10.0^bin/2.0],[Xhist+10.0^bin/2.0]])
;         xfill = reform(xfill, n_elements(xfill))
;         xfill = [xfill[0], xfill, xfill[n_elements(xfill)-1]]
;         yfill = transpose([[yhist],[yhist]])
;         yfill = reform(yfill, n_elements(yfill))
;         yfill = [0, yfill, 0]
          xfill = xfill > 10^!X.CRANGE[0] < 10^!X.CRANGE[1] ;Make sure within plot range
          yfill = yfill > !Y.CRANGE[0] < !Y.CRANGE[1]
       endif else begin
          xfill = transpose([[Xhist-bin/2.0],[Xhist+bin/2.0]])
          xfill = reform(xfill, n_elements(xfill))
          xfill = [xfill[0], xfill, xfill[n_elements(xfill)-1]]
          yfill = transpose([[yhist],[yhist]])
          yfill = reform(yfill, n_elements(yfill))
          yfill = [0, yfill, 0]
          xfill = xfill > !X.CRANGE[0] < !X.CRANGE[1] ;Make sure within plot range
          yfill = yfill > !Y.CRANGE[0] < !Y.CRANGE[1]
       endelse

       if keyword_set(Fcolor) then Fc = im_color(Fcolor) else Fc = !P.Color
       if keyword_set(Fline) then begin
          if keyword_set(Fspacing) then Fs = Fspacing else Fs = 0
          if keyword_set(Forientation) then Fo = Forientation else Fo = 0
          polyfill, xfill,yfill, color=Fc, /line_fill, spacing=Fs, orient=Fo
       endif else begin
          if keyword_set(Fpattern) then begin
             polyfill, xfill,yfill, color=Fc, pattern=Fpattern
          endif else begin
             polyfill, xfill,yfill, color=Fc
          endelse
       endelse

; because the POLYFILL can erase/overwrite parts of the originally plotted
; histogram, we need to replot it here.
       if keyword_set(cumulative) then $
         plot, xhist, yhist, xlog=xlog, _extra=extra  else $
           oplot, [xhist[0]-bin,xhist,xhist[n_hist-1]+bin], [0,yhist,0], $
         psym=psym, _extra=extra
    endif

return    
end
