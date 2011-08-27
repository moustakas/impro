;+
; NAME:
;   IM_MWRFITS
;
; PURPOSE:
;   Simple wrapper on MWRFITS with some bells and whistles.
;
; INPUTS: 
;   outstruct - output structure of image to write out
;   outfile - output file name
;
; OPTIONAL INPUTS: 
;   hdr - output FITS-style header
;
; KEYWORD PARAMETERS: 
;   silent - suppress messages to STDOUT
;   append - append to an existing structure (default is to use
;     /CREATE) 
;   gzip - compress the output file (only needed after using /APPEND) 
;   nogzip - do not compress the output file (default *is* to
;     compress unless /APPEND)
;   clobber - overwrite existing output file, otherwise exit
;
; OUTPUTS: 
;   Writes out a FITS file.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Mar 18, NYU
;   jm09nov10ucsd - added /CLOBBER keyword
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

pro im_mwrfits, outstruct, outfile, hdr, silent=silent, $
  append=append, gzip=gzip, nogzip=nogzip, clobber=clobber

    if keyword_set(append) then begin
       if (file_test(outfile) eq 0) then begin
          splog, 'Output file '+outfile+' does not exist!'
          return
       endif
       if (keyword_set(silent) eq 0) then splog, 'Appending to '+outfile
       mwrfits, outstruct, outfile, hdr, /silent
       if keyword_set(gzip) then spawn, 'gzip -f '+outfile, /sh
    endif else begin
       if im_file_test(repstr(outfile,'.gz','')+'.gz',clobber=clobber) then return
       if (keyword_set(silent) eq 0) then splog, 'Writing '+outfile
       mwrfits, outstruct, outfile, hdr, /create
       if (keyword_set(nogzip) eq 0) then spawn, 'gzip -f '+outfile, /sh
    endelse
    
return
end
    
