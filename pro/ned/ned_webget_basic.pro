;+
; NAME:
;       NED_WEBGET_BASIC
;
; PURPOSE:
;       Collects basic data from NED.
;
; INPUTS:
;       galaxy - NED-compatible list of galaxies [NGALAXY]
;
; OPTIONAL INPUTS:
;       outfile - output file name if WRITE=1
;
; KEYWORD PARAMETERS:
;       write - write the basic data as a binary FITS table 
;
; OUTPUTS:
;       basic - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine is a more convenient drop-in replacement routine
;       for NED's batch-retrieval system.  See also
;       PARSE_NED_BYNAME. 
;
; EXAMPLES:
;       IDL> ned_webget_basic, 'ngc5194', b
;       IDL> struct_print, b
;
; TODO:
;       Parse the list of alternate galaxy names.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Nov 30, NYU - written based on
;          NED_WEBGET_DIAMETERS 
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

pro ned_webget_basic, galaxy1, basic, outfile=outfile, write=write
    
    ngalaxy = n_elements(galaxy1)
    if (ngalaxy eq 0L) then begin
       doc_library, 'ned_webget_basic'
       return
    endif

    if (n_elements(outfile) eq 0L) then outfile = 'ned_webget_basic.fits'
    
; initializez the output structure    

    basic = {$
      galaxy:                       '...', $
      ned_galaxy:                   '...', $
      ra:                           '...', $
      dec:                          '...', $
      z:                           -999.0, $
      z_err:                       -999.0, $
      morph:                        '...', $
      class:                        '...', $
      mag:                          '...', $
;     mag:                         -999.0, $
      dmaj:                        -999.0, $
      dmin:                        -999.0, $
      note:                         '...'}
    basic = replicate(basic,ngalaxy)

    galaxy = strcompress(strupcase(galaxy1),/remove)
    basic.galaxy = galaxy

    for igal = 0L, ngalaxy-1L do begin

       print, format='("Querying NED basic data for galaxy ",I0,"/",I0,".",A10,$)', igal+1L, ngalaxy, string(13b)
       
       nedgalaxy = repstr(galaxy[igal],'+','%2B')

       webinfo = webget('http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname='+$
         nedgalaxy+'&extend=no&search_type=Basic&search_type=Basic')
;      http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=n253&extend=no#BasicData_0
;      webinfo = webget('http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname='+$
;        nedgalaxy+'&extend=no')
       text = webinfo.text
;      openw, lun, 'junk', /get_lun & niceprintf, lun, text & free_lun, lun
       
       noned = where(strmatch(text,'*not currently recognized*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif

       noned = where(strmatch(text,'*there is no object with this name*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif

       noned = where(strmatch(text,'*object name required*',/fold),nnoned)
       if (nnoned ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' not recognized by NED.'
          continue
       endif

       mult = where(strmatch(text,'*multiple choices*',/fold),nmult)
       if (nmult ne 0L) then begin
          splog, 'Galaxy '+galaxy[igal]+' contains multiple choices in NED:'
          istart = where(strmatch(text,'*<PRE>*',/fold) eq 1B,nstart)+1L
          iend = where(strmatch(text,'*</PRE>*',/fold) eq 1B,nend)-1L
          if (nstart ge 1L) and (nend ge 1L) then niceprint, strtrim(text[istart[0]:iend[0]],2)
          continue
       endif

; crop the text
       istart = where(strmatch(text,'*searching ned for object*',/fold) eq 1B,nstart)
       iend = where(strmatch(text,'*back to the list*',/fold) eq 1B,nend)
       if (nstart ge 1L) and (nend ge 1L) then text = text[istart[0]:iend[0]]

; NED name and coordinates
       istart = where(strmatch(text,'*objects found in ned*',/fold) eq 1B,nstart)
       iend = where(strmatch(text,'*detailed information for each object*',/fold) eq 1B,nend)
       if (nstart ne nend) then message, 'Problem, Dave.'

       subtext = text[istart:iend]
       line = (where(strmatch(subtext,'*ObjNo1*',/fold) eq 1B,nline))[0]
       if (nline eq 0L) then message, 'This should not happen!!'

       basic[igal].ned_galaxy = strcompress(strmid(subtext[line],34,30),/remove)

; the NED HTML output is sometimes flaky, and does not return the
; expected code/format, which screws up the galaxy name and the
; coordinates; so create a special case
       if strmatch(basic[igal].ned_galaxy,'*HREF*',/fold) then begin
          galaxy_line = strmid(subtext[line+1L],strpos(subtext[line+1L],'>')+1L,100)
          ra_line = strmid(subtext[line+2L],strpos(subtext[line+2L],'>')+1L,100)
          dec_line = strmid(subtext[line+3L],strpos(subtext[line+3L],'>')+1L,100)
          basic[igal].ned_galaxy = strcompress(galaxy_line,/remove) ; not very robust!!
          basic[igal].ra = repstr(repstr(repstr(strcompress(ra_line,/remove),'h',':'),'m',':'),'s','')
          basic[igal].dec = repstr(repstr(repstr(strcompress(dec_line,/remove),'d',':'),'m',':'),'s','')
       endif else begin
          basic[igal].ra = repstr(repstr(repstr(strcompress(strmid(subtext[line],67,12),$
            /remove),'h',':'),'m',':'),'s','')
          basic[igal].dec = repstr(repstr(repstr(strcompress(strmid(subtext[line],78,12),$
            /remove),'d',':'),'m',':'),'s','')
       endelse

; check if there is an "essential note" for the galaxy       
       yesnote = where(strmatch(basic[igal].ned_galaxy,"*"),nyesnote)
       if nyesnote then begin
          basic[igal].ned_galaxy = repstr(basic[igal].ned_galaxy,'*','')
          line = where((strmatch(text,'*ESSENTIAL NOTE*') eq 1B) and $
            (strmatch(text,'*=>*',/fold) eq 0B),nline)
          if (nline ne 0L) then begin
             if (strmatch(text[line+2],'*N/A*') eq 0) then $
               basic[igal].note = strtrim(text[line+2L],2)
          endif
       endif
       
; parse the basic data        
       istart = where(strmatch(text,'*BASIC DATA*') eq 1B,nstart)
       iend = where(strmatch(text,'*QUANTITIES*') eq 1B,nend)
       if (nend eq 0L) then iend = where(strmatch(text,'*back to the list*',/fold) eq 1B,nend)
       if (nstart ne nend) then message, 'Problem, Dave.'

       subtext = text[istart:iend-1L]

; redshift and error
       line = where(strmatch(subtext,'*Redshift*:*') eq 1B,nline)
       if (nline ne 0L) then begin
          i1 = strpos(subtext[line],':')+1
          i2 = strpos(subtext[line],'+/-')
          if (i2[0] ne -1L) then begin
             basic[igal].z = strmid(subtext[line],i1,i2-i1)
             basic[igal].z_err = strmid(subtext[line],i2+3,15)
             if (basic[igal].z_err le 0.0) then basic[igal].z_err = -999.0
          endif else begin
             basic[igal].z = strmid(subtext[line],i1,20)
          endelse
       endif

; diameters
       line = where(strmatch(subtext,'*Major Diameter*:*') eq 1B,nline)
       if (nline ne 0L) then begin
          dmaj = strmid(subtext[line],strpos(subtext[line],':')+1,30)
          if (strcompress(dmaj,/remove) ne '') then basic[igal].dmaj = dmaj
       endif

       line = where(strmatch(subtext,'*Minor Diameter*:*') eq 1B,nline)
       if (nline ne 0L) then begin
          dmin = strmid(subtext[line],strpos(subtext[line],':')+1,30)
          if (strcompress(dmin,/remove) ne '') then basic[igal].dmin = dmin
       endif

; magnitude
       line = where(strmatch(subtext,'*Magnitude*:*') eq 1B,nline)
       if (nline ne 0L) then begin
          mag = strmid(subtext[line],strpos(subtext[line],':')+1,30)
          if (strcompress(mag,/remove) ne '') then $
            basic[igal].mag = mag
       endif

; morphology & classifications; if there is a semi-colon then assume
; that NED provides a spectral classification
       line = where(strmatch(subtext,'*Classifications*:*') eq 1B,nline)
       if (nline ne 0L) then begin

          i1 = strpos(subtext[line],':')+1

          if strmatch(subtext[line],'*;*') then begin
             i2 = strpos(subtext[line],';')
             basic[igal].class = strcompress(strmid(subtext[line],i2+1,30),/remove)
          endif else begin
             i2 = i1+30
          endelse

          morph = strcompress(strmid(subtext[line],i1,i2-i1),/remove)
          if (strcompress(morph,/remove) ne '') then $
            basic[igal].morph = morph

          if strmatch(basic[igal].morph,'*LINER*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'LINER','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'LINER' else $
               basic[igal].class = basic[igal].class+'LINER'
          endif
          if strmatch(basic[igal].morph,'*NLAGN*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'NLAGN','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'NLAGN' else $
               basic[igal].class = basic[igal].class+'NLAGN'
          endif
          if strmatch(basic[igal].morph,'*Sy1*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'Sy1','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'Sy1' else $
               basic[igal].class = basic[igal].class+'Sy1'
          endif
          if strmatch(basic[igal].morph,'*Sy2*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'Sy2','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'Sy2' else $
               basic[igal].class = basic[igal].class+'Sy2'
          endif
          if strmatch(basic[igal].morph,'*Sy1.9*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'Sy1.9','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'Sy1.9' else $
               basic[igal].class = basic[igal].class+'Sy1.9'
          endif
          if strmatch(basic[igal].morph,'*Sy*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'Sy','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'Sy' else $
               basic[igal].class = basic[igal].class+'Sy'
          endif
          if strmatch(basic[igal].morph,'*HII*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'HII','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'HII' else $
               basic[igal].class = basic[igal].class+'HII'
          endif
          if strmatch(basic[igal].morph,'*Sbrst*') then begin
             basic[igal].morph = repstr(basic[igal].morph,'Sbrst','')
             if (strtrim(basic[igal].class,2) eq '...') then basic[igal].class = 'Sbrst' else $
               basic[igal].class = basic[igal].class+'Sbrst'
          endif

       endif
          
       delvarx, webinfo, text
;      help, basic[igal], /str & stop
       
    endfor

; write out    
    
    if keyword_set(write) then begin
       splog, 'Writing '+outfile
       mwrfits, basic, outfile, /create
       spawn, 'gzip -f '+outfile, /sh
    endif
    
return
end
    
