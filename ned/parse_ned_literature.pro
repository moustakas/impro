;+
; NAME:
;   PARSE_NED_LITERATURE
;
; PURPOSE:
;   Parse the output of a NED literature search.
;
; INPUTS: 
;   nedfile - literature results file from NED
;
; OPTIONAL INPUTS: 
;   outfile - output file name (default 'outfile.fits' or
;     'outfile.txt'; see /TEXTFILE)
;
; KEYWORD PARAMETERS: 
;   textfile - write out an ASCII file (default is to write a binary
;     FITS table)
;   kms - convert redshifts to km/s (fragile!)
;
; OUTPUTS: 
;   Writes out a FITS or ASCII file with the derived literature
;   results. 
;
; COMMENTS:
;   Sorts the output by RA.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Nov 07, U of A
;   jm02mar22uofa - modified to output a binary fits table; any
;     existing binary fits table of the same name is overwritten 
;
; Copyright (C) 2001-2002, John Moustakas
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

pro parse_ned_literature, nedfile, outfile=outfile, textfile=textfile, kms=kms

    if n_elements(nedfile) eq 0 then begin
       doc_library, 'parse_ned_literature'
       return
    endif
    
; find out how many objects are in the file
    if file_test(nedfile,/regular) eq 0L then begin
       splog, 'File '+nedfile+' not found.'
       return
    endif
    
    objindx = 19
    nobjline = djs_readilines(nedfile,indx=objindx)
    
    reads, nobjline, nobj
    nobj = long(nobj)
    splog, 'There are '+strn(nobj)+' objects in '+nedfile+'.'
    
    lstart = 22L                ; starting line
    lend = lstart + 2*nobj-1L   ; ending line
    lindx = lindgen((lend-lstart)+1)+lstart
    
    nedata = djs_readilines(nedfile,indx=lindx)
    
    ginfo = strarr(8,nobj)
    data = create_struct('galaxy', '', 'ra', '', 'dec', '', 'z', -999D, 'z_err', -999D)
    data = replicate(data,nobj)

    if keyword_set(kms) then zscale = 2.99792458E5 else zscale = 1.0
    
    for i = 0L, nobj-1L do begin
       
       line = strjoin(nedata[i*2:i*2+1])

       data[i].galaxy = strcompress(strmid(line,7,16)+strmid(line,89,12),/remove) ; galaxy

       data[i].ra = strjoin(strsplit(strmid(line,24,11),'hms',/extr),':')
       data[i].dec = strjoin(strsplit(strmid(line,36,10),'dms',/extr),':')
       
       z = strcompress(strmid(line,46,9),/remove)
       zerr = strcompress(strmid(line,126,9),/remove)

       if (z ne '') then data[i].z = double(z)/zscale
       if (zerr ne '') then if (float(zerr) ne 0.0) then data[i].z_err = double(zerr)/zscale
       
    endfor
    
; sort by RA
    srtra = sort(im_hms2dec(data.ra))
    data = data[srtra]

; write a binary fits table or textfile
    if keyword_set(textfile) then begin
       if not keyword_set(outfile) then outfile = 'outfile.txt'
       splog, 'Writing '+outfile+'.'
       print_struct, data, file=outfile
    endif else begin
       if not keyword_set(outfile) then outfile = 'outfile.fits'
       splog, 'Writing '+outfile+'.'
       mwrfits, data, outfile, /create
       spawn, ['gzip -f '+outfile], /sh
    endelse
    
return
end
