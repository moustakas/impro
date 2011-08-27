;+
; NAME:
;      PARSE_NED_BYNAME
;
; PURPOSE:
;       Parse basic galaxy properties from a NED batch request.
;
; INPUTS:
;       nedfile      - basic data text file from NED
;
; OPTIONAL INPUTS:
;       inputnedfile - input text file sent to NED with a list of
;                      objects, one per line (spaces and funny
;                      characters are okay)
;       nedpath      - path name to NEDFILE and INPUTNEDFILE
;       outfile      - output file name (either OUTFILE.FITS or
;                      OUTFILE.TXT, depending on whether TEXTFILE=1)
;       outpath      - output path name for OUTFILE
;
; KEYWORD PARAMETERS:
;       sortra       - sort the output data structure by RA
;       textfile     - write a columnular text file rather than a
;                      binary FITS file
;
; OUTPUTS:
;       data         - output photometry data structure
;          GALAXY        - input galaxy name from INPUTNEDFILE
;          NEDGALAXY     - NED galaxy name from NEDFILE
;          RA            - right ascenscion (J2000)
;          DEC           - declination (J2000)
;          Z             - redshift (dimensionless)
;          Z_ERR         - redshift error (dimensionless)
;          MORPH         - NED morphology
;          MAG           - "optical" magnitude (not reliable)
;          DMAJ          - major axis length at 25 mag/arcsec2 (arcmin)
;          DMIN          - minor axis length at 25 mag/arcsec2 (arcmin)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;      Any existing binary fits table of the same name is
;      overwritten.  
;
; EXAMPLE:
;       IDL> parse_ned_byname, 'ned.dat', outfile='ned.fits'
;       IDL> parse_ned_byname, 'ned.dat', /textfile, /sortra, outfile='ned.txt'
;
; MODIFICATION HISTORY:
;      J. Moustakas, 2001 November 7, U of A
;      jm03apr3uofa, cleaned up and documented
;      jm03sep17uofa, added error checking if the object does not
;         exist in NED
;      jm06jan24uofa - added UNIQUE keyword
;      jm07nov28nyu - check whether the NED redshift is dimensionless
;                     or in km/s, and convert as appropriate
;
; Copyright (C) 2001, 2003, 2006-2007, John Moustakas
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

pro parse_ned_byname, nedfile, data, inputnedfile=inputnedfile, nedpath=nedpath, $
  sortra=sortra, unique=unique, textfile=textfile, outfile=outfile, outpath=outpath
    
    if (n_elements(nedfile) eq 0L) then begin
       doc_library, 'parse_ned_byname'
       return
    endif
    
    cspeed = 299850D ; this is what NED uses [km/s]
;   cspeed = 2.99792458D5 ; [km/s]

    if (n_elements(nedpath) eq 0L) then nedpath = cwd()
    
; read the whole file and crop

    if file_test(nedpath+nedfile,/regular) eq 0L then begin
       splog, 'File '+nedpath+nedfile+' not found.'
       return
    endif
    
    nedata = djs_readilines(nedpath+nedfile)
    lstart = (where(strmatch(nedata,'*SEARCH RESULTS*') eq 1B))[0]
    lend = (where(strmatch(nedata,'*all the objects found*') eq 1B))[0]

    nedata = nedata[lstart:lend]
    
; now find out how many objects are in the file by doing a search on BYNAME

    objname = where(strmatch(nedata,'BYNAME *') eq 1B,nobj)
    splog, 'There are '+strn(nobj)+' objects in '+nedfile+'.'

; read the input text file to NED, if it was given

    if n_elements(inputnedfile) ne 0L then begin

       if file_test(nedpath+inputnedfile,/regular) eq 0L then begin
          splog, 'Input ned file '+nedpath+inputnedfile+' not found.'
          return
       endif

       inputfile = strcompress(djs_readlines(nedpath+inputnedfile),/remove)
       galaxy = inputfile[where((inputfile ne '') and (strmatch(inputfile,'##*') eq 0B))]
       ngalaxy = n_elements(galaxy)

       splog, 'There are '+strn(ngalaxy)+' objects in '+inputnedfile+'.'

       if (nobj ne ngalaxy) then begin
          splog, 'NED did not return basic data for all the objects . . . returning.'
          return
       endif

    endif
    
    data = {$
      galaxy:            ' ', $
      nedgalaxy:         ' ', $
      ra:       '99:99:99.0', $
      dec:       '+99:99:99', $
      z:               -999D, $
      z_err:           -999D, $
      morph:             ' ', $
      mag:             -999.0,$
      dmaj:            -999.0,$
      dmin:            -999.0}
    data = replicate(data,nobj)

    for i = 0L, nobj-1L do begin

; check to see if none or one object was found

       searchgalaxy = strcompress(strmid(nedata[objname[i]],6,50),/remove)
       if strmatch(nedata[objname[i]+1],'*does not exist*',/fold) then begin
          splog, 'Object '+searchgalaxy+' does not exist in the data base.'
          data[i].galaxy = searchgalaxy
       endif
       
       if strmatch(nedata[objname[i]+1],'*1*object*found*',/fold) then begin

          line = strjoin(nedata[objname[i]+4:objname[i]+5])

          if n_elements(ngalaxy) ne 0L then data[i].galaxy = strcompress(galaxy[i],/remove) else $
            data[i].galaxy = searchgalaxy
          data[i].nedgalaxy = strcompress(strmid(line,7,16)+strmid(line,89,12),/remove) ; NED galaxy name

          data[i].ra = strjoin(strsplit(strmid(line,24,11),'hms',/extr),':')
          data[i].dec = strjoin(strsplit(strmid(line,36,10),'dms',/extr),':')

          z = strcompress(strmid(line,46,9),/remove)
          zerr = strcompress(strmid(line,126,9),/remove)

          if (z ne '') then data[i].z = z
          if (zerr ne '') then data[i].z_err = zerr

          data[i].morph  = strcompress(strmid(line,60,20),/remove) ; morphology

          mag = strcompress(strmid(line,140,8),/remove)
          if (mag ne '') and (mag ne '/?') then data[i].mag = float(mag)
          
          dmaj = strcompress(strmid(line,148,5),/remove)
          if dmaj ne '' then data[i].dmaj = float(dmaj)

          dmin = strcompress(strmid(line,155,5),/remove)
          if dmin ne '' then data[i].dmin = float(dmin)

       endif
          
    endfor

; check and see if the redshift was returned in km/s or as "z"; the
; way I do this check is pretty dumb

    cz = where(data.z gt 10.0,ncz)
    if (ncz ne 0L) then begin
       good = where((data.z gt -900.0),ngood)
       if (ngood ne 0L) then data[good].z = data[good].z/cspeed
       good = where((data.z_err gt -900.0),ngood)
       if (ngood ne 0L) then data[good].z_err = data[good].z_err/cspeed
    endif
    
; sort by RA

    if keyword_set(sortra) then begin
       srtra = sort(im_hms2dec(data.ra))
       data = data[srtra]
    endif

; only include unique objects

    if keyword_set(unique) then begin
       srt = uniq(strtrim(data.nedgalaxy,2),sort(strtrim(data.nedgalaxy,2)))
       data = data[srt]
       splog, 'Writing '+string(n_elements(data),format='(I0)')+' unique objects.'
    endif

; write a binary fits table
    
    if not keyword_set(outpath) then outpath = cwd()
    if keyword_set(textfile) then begin
       if not keyword_set(outfile) then outfile = 'outfile.txt'
       splog, 'Writing '+outpath+outfile+'.'
       print_struct, data, file=outpath+outfile
    endif else begin
       if not keyword_set(outfile) then outfile = 'outfile.fits'
       splog, 'Writing '+outpath+outfile+'.gz.'
       mwrfits, data, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh
    endelse
    
return
end
