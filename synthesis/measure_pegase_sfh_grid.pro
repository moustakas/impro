;+
; NAME:
;       MEASURE_PEGASE_SFH_GRID
;
; PURPOSE:
;       Measure quantities of interest in the Pegase SFH grid and
;       convert to FITS files.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       gridpath - output path [default getenv('PEGASE_HR_SFHGRID_DIR')] 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Should be run immediately after WRITE_PEGASE_SFH_GRID. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 19, U of A - based on WRITE_SFH_BASE_MODELS 
;       jm07feb26nyu - update to used PEGASE-HR
;
; Copyright (C) 2006-2007, John Moustakas
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

pro measure_pegase_sfh_grid, data, write=write

    gridpath = getenv('PEGASE_HR_SFHGRID_DIR')+'/' ; SFH grid path

; use all the ages - jm08apr20nyu    
    
;   agemax = ceil(getage(0.0))  ; ceil(getage(0.0)*10.0)/10.0
;   agemin = 0.01 & dage = 0.01
;   age = findgen((agemax-agemin)/dage+1)*dage+agemin
;   age = age*1E3
    
;   peglist = file_basename(file_search(gridpath+'*tau_000.0*.fits',count=pegcount))
    peglist = file_basename(file_search(gridpath+'*.fits',count=pegcount))
    if (pegcount eq 0L) then begin
       splog, 'No PEGASE files found in '+gridpath
       return
    endif

    fitslist = repstr(peglist,'.fits','.info.fits')
    
    for ii = 0L, pegcount-1L do begin
       peg = im_read_peg(gridpath+peglist[ii])
       data = im_measure_peg(peg,age=age)
       if keyword_set(write) then begin
          splog, 'Writing '+gridpath+'measure/'+fitslist[ii]
          mwrfits, data, gridpath+'measure/'+fitslist[ii], /create
       endif
    endfor

return
end
    
