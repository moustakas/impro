;+
; NAME:
;   IM_HEADERFORAGE()
;
; PURPOSE:
;   Return the contents of a FITS file header as a structure.
;
; INPUTS: 
;   fitslist - list of FITS files
;
; OPTIONAL INPUTS: 
;   ext - extension number to read
;
; OUTPUTS: 
;   forage - data structure with all the header entries as individual
;     structure tags
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Apr 07, U of A
;   jm08aug06nyu - removed DATAPATH input; absolute file name now
;     assumed; EXT optional input added 
;
; Copyright (C) 2005, 2008, John Moustakas
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

function im_headerforage, fitslist, ext=ext

    if (n_elements(fitslist) eq 0L) then begin
       print, 'Syntax - forage = im_headerforage()'
       return, -1
    endif

    fitslist = file_search(fitslist,count=nfits)
    if (n_elements(ext) eq 0L) then ext = 0L
    next = n_elements(ext)

    for jj = 0L, nfits-1L do begin

       for ii = 0L, next-1L do begin
          
          h = headfits(fitslist[jj],ext=ext[ii])
          s = im_hdr2struct(h)
          forage1 = struct_trimtags(s,except='*HISTORY*')

; add tags here

          add = {fitsfile: fitslist[jj], extension: ext[ii]}
          forage1 = struct_addtags(add,forage1)

          if (n_elements(forage) eq 0L) then forage = forage1 else begin
             forage = struct_trimtags(forage,select=tag_names(forage1))
             forage1 = struct_trimtags(forage1,select=tag_names(forage[0]))
             forage = [forage,forage1]
          endelse

       endfor
          
    endfor

return, reform(forage)
end    
