;+
; NAME:
;   GALEX_RESOLVE_CASJOBS_DUPLICATES
;
; PURPOSE:
;   Analyze the output of a CasJobs query to resolve/remove
;   duplicates. 
;
; INPUTS: 
;   incat - input catalog casjobs catalog
;
; OPTIONAL INPUTS: 
;   idtagname - tag name to use to identify repeat observations of the
;     same object (see EXAMPLES, below)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   IDL> in = mrdfits(ages_path(/mycat)+'galex/ages_galex_gr45_casjobs.fits.gz'
;   IDL> out = galex_resolve_casjobs_duplicates(in,idtag='ages_id')
;   IDL> best = out[where(out.isbest)]
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Apr 30, UCSD
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

function galex_resolve_casjobs_duplicates, incat, idtagname=idtagname

    nobj = n_elements(incat)
    if (nobj eq 0L) then begin
       doc_library, 'galex_resolve_casjobs_duplicates'
       return, -1
    endif

    if (n_elements(idtagname) eq 0) then idtagname = 'id'
    idtag = tag_indx(incat,idtagname)
    if (idtag[0] eq -1) then begin
       splog, 'ID tag name '+idtagname+ 'not found!'
       return, -1
    endif
    
    addtags = replicate({ndup: 0, isbest: 1},nobj)
    outcat = struct_addtags(incat,addtags)

    allid = incat.(idtag)
    id = allid[uniq(allid,sort(allid))]
    nid = n_elements(id)
    if (nid eq n_elements(allid)) then begin
       splog, 'No duplicate observations - woo hoo!'
       return, outcat
    endif
    
    for ii = 0L, nid-1L do begin
       print, "Working on object ", ii, " of ", nid, string(13b), $
         format='(A0,I0,A0,I0,A1,$)'
       these = where(id[ii] eq allid,nthese)
; check for repeat observations; take the one with the highest S/N,
; first in the NUV channel, and then in the FUV channel
       if (nthese gt 1L) then begin 
          outcat[these].ndup = nthese
          outcat[these].isbest = 0
; check NUV first
          nuv_good = where(outcat[these].nuv_fluxerr_auto gt 0.0,nnuv_good)
          if (nnuv_good ne 0) then begin
             nuv_snr = outcat[these[nuv_good]].nuv_flux_auto/outcat[these[nuv_good]].nuv_fluxerr_auto
;            nuv_snr = 2.5/(alog(10)*outcat[these[nuv_good]].nuv_magerr)
             nuv_snrmax = where(nuv_snr eq max(nuv_snr),nnuvsnrmax)
             if (nnuvsnrmax eq 1) then outcat[these[nuv_good[nuv_snrmax]]].isbest = 1 else begin
; sometimes the magnitude error (i.e., the S/N) is identical!?!
                splog, 'Identical NUV S/N for object ID '+strtrim(id[ii],2)+'!'
                outcat[these[nuv_good[nuv_snrmax[0]]]].isbest = 1 ; pick the first one
             endelse
          endif else begin ; check FUV
             fuv_good = where(outcat[these].fuv_fluxerr_auto gt 0.0,nfuv_good)
;            fuv_good = where(outcat[these].fuv_mag gt 0.0,nfuv_good)
             if (nfuv_good ne 0) then begin
                fuv_snr = outcat[these[fuv_good]].fuv_flux_auto/outcat[these[fuv_good]].fuv_fluxerr_auto
;               fuv_snr = 1.0/(alog(10)*outcat[these[fuv_good]].fuv_magerr)
                fuv_snrmax = where(fuv_snr eq max(fuv_snr),nfuvsnrmax)
                if (nfuvsnrmax eq 1) then outcat[these[fuv_good[fuv_snrmax]]].isbest = 1 else begin
                   splog, 'Identical FUV S/N for object ID '+strtrim(id[ii],2)+'!'
                   outcat[these[fuv_good[fuv_snrmax[0]]]].isbest = 1 ; pick the first one
                endelse
             endif
          endelse
          check = where(outcat[these].isbest eq 1)
          if (check[0] eq -1) then message, 'Fix me!'
       endif 
    endfor

return, outcat
end
