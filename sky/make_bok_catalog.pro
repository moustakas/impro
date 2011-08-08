;+
; NAME:
;   MAKE_BOK_CATALOG
;
; PURPOSE:
;   Generate a catalog that for the Bok/90-inch.
;
; INPUTS: 
;   ra, dec - input celestial coordinates
;   mag - magnitude
;   object - object name
;   epoch - coordinate epoch
;
; OPTIONAL INPUTS: 
;   catname - output catalog name (default 'catalog')
;
; OUTPUTS: 
;   Writes out an appropriately formated catalog.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Jun 24, U of A
;
; Copyright (C) 2002, John Moustakas
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

pro make_bok_catalog, ra, de, mag, object, epoch, catname=catname

    nobject = n_elements(ra)
    id = lindgen(nobject)
    
    if n_elements(catname) eq 0L then catname = 'catalog'
    
    dtor = !dpi/180.0D ; degrees --> radians
    
; convert RA and DEC to radians

    rra = 15D0*dtor*im_hms2dec(ra) ; [rad]
    rde = dtor*im_hms2dec(de)      ; [rad]

; define the output formats

    frmtpos = '(x,I3,x,I4,x,F12.10,x,A1,F12.10,6x,A1,6x,A1,x,A9,56x,F7.2)'
    frmtneg = '(x,I3,x,I4,x,F12.10,x,F13.10,6x,A1,6x,A1,x,A9,56x,F7.2)'

    rdepos = where(rde gt 0.0,nrdepos,comp=rdeneg,ncomp=nrdeneg)
    
    openw, lun, catname+'.cat', /get_lun
    for i = 0L, nobject-1L do begin
       if rde[i] gt 0.0 then printf, lun, i, mag[i], rra[i], '+', $
         rde[i], '0', '0', object[i], epoch[i], format=frmtpos else $
         printf, lun, i, mag[i], rra[i], rde[i], '0', '0', $
         object[i], epoch[i], format=frmtneg
    endfor
    free_lun, lun

return
end
