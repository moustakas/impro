;+
; NAME:
;   IM_READ_STARBURST99()
;
; PURPOSE:
;   Read the Starburst99 spectra and data files into a structure. 
;
; INPUTS:
;   specfile - name of the *spectrum* output file
;
; OPTIONAL INPUTS:
;   ewidthfile - name of the *ewidth* output file
;   quantafile - name of the *quanta* output file
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   sb99 - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   See http://www.stsci.edu/science/starburst99/docs/run.html for a
;   description of the various output files.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 12, UCSD
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

function im_read_starburst99, specfile, ewidthfile=ewidthfile, $
  quantafile=quantafile

    if (n_elements(specfile) eq 0) then begin
       doc_library, 'im_read_starburst99'
       return, -1
    endif

    if (file_test(specfile) eq 0) then begin
       print, 'Spectrum file '+specfile+' not found'
       return, -1
    endif

; read the spectrum file, parse, and pack into a data structure
    splog, 'Reading '+specfile
    readfast, specfile, data, skip=7, /double
    alltime = reform(data[0,*])
    time = alltime[uniq(alltime,sort(alltime))]
    ntime = n_elements(time)
    for ii = 0L, ntime-1L do begin
       these = where(time[ii] eq alltime,npix)
       if (ii eq 0) then begin
          sb99 = {$
            time: 0.0D,$
            wave:      dblarr(npix),$
            lum_tot:   dblarr(npix),$
            lum_stars: dblarr(npix),$
            lum_gas:   dblarr(npix)}
          sb99 = replicate(sb99,ntime)
       endif
       sb99[ii].time = time[ii]
       sb99[ii].wave = reform(data[1,these])
       sb99[ii].lum_tot = reform(data[2,these])
       sb99[ii].lum_stars = reform(data[3,these])
       sb99[ii].lum_gas = reform(data[4,these])
    endfor

; check for the *ewidth* file
    if (n_elements(ewidthfile) ne 0) then begin
       if (file_test(ewidthfile) eq 0) then begin
          print, 'Spectrum file '+ewidthfile+' not found'
       endif else begin
          moretags = {$
            ha_c:  0.0, ha:  0.0, ewha:  0.0,$ ; H-alpha
            hb_c:  0.0, hb:  0.0, ewhb:  0.0,$ ; H-beta
            pab_c: 0.0, pab: 0.0, ewpab: 0.0,$ ; Pa-beta
            brg_c: 0.0, brg: 0.0, ewbrg: 0.0}  ; Br-gamma
          sb99 = struct_addtags(temporary(sb99),replicate(moretags,ntime))
          splog, 'Reading '+ewidthfile
          readfast, ewidthfile, data, skip=7
          newtime = reform(data[0,*]) ; this file has much higher time resolution
; H-alpha
          sb99.ha_c = interpol(reform(data[1,*]),newtime,sb99.time)
          sb99.ha   = interpol(reform(data[2,*]),newtime,sb99.time)
          sb99.ewha = interpol(reform(data[3,*]),newtime,sb99.time)
; H-beta
          sb99.hb_c = interpol(reform(data[4,*]),newtime,sb99.time)
          sb99.hb   = interpol(reform(data[5,*]),newtime,sb99.time)
          sb99.ewhb = interpol(reform(data[6,*]),newtime,sb99.time)
; Pa-beta
          sb99.pab_c = interpol(reform(data[7,*]),newtime,sb99.time)
          sb99.pab   = interpol(reform(data[8,*]),newtime,sb99.time)
          sb99.ewpab = interpol(reform(data[9,*]),newtime,sb99.time)
; Br-gamma
          sb99.brg_c = interpol(reform(data[10,*]),newtime,sb99.time)
          sb99.brg   = interpol(reform(data[11,*]),newtime,sb99.time)
          sb99.ewbrg = interpol(reform(data[12,*]),newtime,sb99.time)
       endelse
    endif

; check for the *quanta* file (not fully parsed!)
    if (n_elements(quantafile) ne 0) then begin
       if (file_test(quantafile) eq 0) then begin
          print, 'Spectrum file '+quantafile+' not found'
       endif else begin
          moretags = {nhi: 0.0, nhei: 0.0, nheii: 0.0}
          sb99 = struct_addtags(temporary(sb99),replicate(moretags,ntime))
          splog, 'Reading '+quantafile
          readfast, quantafile, data, skip=7
          newtime = reform(data[0,*]) ; this file has much higher time resolution
          sb99.nhi   = interpol(reform(data[1,*]),newtime,sb99.time)
          sb99.nhei  = interpol(reform(data[3,*]),newtime,sb99.time)
          sb99.nheii = interpol(reform(data[5,*]),newtime,sb99.time)
       endelse
    endif
          
return, sb99
end
