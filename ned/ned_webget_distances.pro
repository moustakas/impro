;+
; NAME:
;   NED_WEBGET_DISTANCES
;
; PURPOSE:
;   Collect distance references from NED.
;
; INPUTS:
;   gal - NED-compatible list of galaxies [NGALAXY]
;
; OPTIONAL INPUTS:
;   outfile - output file name if WRITE=1
;
; KEYWORD PARAMETERS:
;   write - write DIST as a binary FITS table
;
; OUTPUTS:
;   dist - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Sep 22, UCSD
;
; Copyright (C) 2009, John Moustakas
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

pro ned_webget_distances, gal, dist, outfile=outfile, write=write
    
    ngal = n_elements(gal)
    if (ngal eq 0L) then begin
       doc_library, 'ned_webget_dist'
       return
    endif

    if (n_elements(outfile) eq 0) then outfile = 'ned_webget_dist.fits'
    
; initializez the output structure    
    delvarx, dist
    dist_template = {$
      galaxy:         '', $
;     ned_galaxy:     '', $
      dmod:       -999.0, $
      dmod_err:   -999.0, $
      distance:   -999.0, $
      method:         '', $
      refcode:        '', $
      notes:          ''}
    
    for igal = 0L, ngal-1L do begin

       nedgal = repstr(gal[igal],'+','%2B')
       webinfo = webget('http://nedwww.ipac.caltech.edu/cgi-bin/nDistance?name='+nedgal)
       text = webinfo.text

       dist1 = dist_template
       dist1.galaxy = strcompress(strupcase(gal[igal]),/remove)

       line = where(strmatch(text,'*distances found in ned*',/fold),nline)
       ndist = fix(strmid(text[line],4,2))
       if (ndist eq 0) then begin
          splog, 'No distances found for object '+gal[igal]
          if (n_elements(dist) eq 0) then dist = dist1 else $
            dist = [dist,dist1]
          continue
       endif

       dist1 = replicate(dist1,ndist)
       
       istart = where(strmatch(text,'*individually referenced*',/fold) eq 1B,nstart)
       iend = where(strmatch(text,'*back to ned*',/fold) eq 1B,nend)
       if (nstart ne nend) then message, 'Problem, Dave.'

       subtext = text[istart+1:iend-1]

       for kk = 0, n_elements(subtext)-1 do begin
          if strmatch(subtext[kk],'*<tr><td>*') and $
            strmatch(subtext[kk],'*</td></tr>*') then begin
             i1 = strpos(subtext[kk],'<tr><td>')
             i2 = strpos(subtext[kk],'</td></tr>')
             info = strmid(subtext[kk],i1,i2-i1)
             
          endif
       endfor
       
       help, subtext

    endfor
       
stop       
; write out    
    if keyword_set(write) then begin
       splog, 'Writing '+outfile
       mwrfits, dist, outfile, /create
       spawn, ['gzip -f '+outfile], /sh
    endif
    
return
end
    
