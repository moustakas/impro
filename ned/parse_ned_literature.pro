;+
; NAME:
;   PARSE_NED_LITERATURE
;
; PURPOSE:
;   Parse 
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
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
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



pro parse_ned_literature, nedfile, outfile=outfile, outpath=outpath, $
  textfile=textfile, kms=kms
; jm01nov7uofa
; this routine will parse a ned text file from a literature search
; jm02mar22uofa - modified to output a binary fits table; any existing
;                 binary fits table of the same name is overwritten

; if the redshifts are in km/s then set this keyword - doesn't work perfectly!!
    
    if n_elements(nedfile) eq 0L then begin
       print, 'Syntax - parse_ned_literature, nedfile, [outfile=, outpath=], textfile=textfile, kms=kms'
       return
    endif
    
    cspeed = 2.99792458E5 ; [km/s]

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

    if keyword_set(kms) then zscale = 2.99D5 else zscale = 1.0
    
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
    
    if not keyword_set(outpath) then outpath = cwd()
    if keyword_set(textfile) then begin
       if not keyword_set(outfile) then outfile = 'outfile.txt'
       splog, 'Writing '+outpath+outfile+'.'
       print_struct, data, file=outpath+outfile
    endif else begin
       if not keyword_set(outfile) then outfile = 'outfile.fits'
       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, data, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh
    endelse
    
return
end
