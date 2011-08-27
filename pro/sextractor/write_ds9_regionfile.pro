;+
; NAME:
;   WRITE_DS9_REGIONFILE
;
; PURPOSE:
;   Write out a DS9 region file.
;
; INPUTS: 
;   ra  - right ascension
;   dec - declination 
;
; OPTIONAL INPUTS: 
;   filename - output file name (default 'ds9.reg')
;   symbol   - symbol to draw; supported symbols are: 'box',
;              'diamond', 'cross', 'x', 'arrow', and 'boxcircle'
;              (default 'circle') 
;   color    - region color; supported colors are 'white',
;              'black', 'red', 'green', 'blue', 'cyan', 'magenta',
;              and 'yellow' (default 'green')
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
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 18, NYU - written
;   jm08jun09nyu - added SYMBOL optional input
;
; Copyright (C) 2007-2008, John Moustakas
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

pro write_ds9_regionfile, ra, dec, filename=filename, symbol=symbol, $
  color=color

    nobj = n_elements(ra)
    if (nobj eq 0L) or (nobj ne n_elements(dec)) then begin
       doc_library, 'write_ds9_regionfile'
       return
    endif
    
    if (n_elements(filename) eq 0L) then filename = 'ds9.reg'
    if (n_elements(symbol) eq 0L) then symbol = 'circle'
    if (n_elements(color) eq 0L) then color = 'green'
    
    openw, lun, filename, /get_lun
    printf, lun , 'global color='+color+' font="helvetica 10 normal" '+$
      'select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
    printf, lun, 'fk5'
    for ii = 0L, nobj-1L do printf, lun, 'point('+string(ra[ii],$
      format='(E16.10)')+','+string(dec[ii],format='(E17.10)')+$
      ') # point='+strtrim(symbol,2)
    free_lun, lun

return
end
    
