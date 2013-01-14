;+
; NAME:
;   READ_SFHGRID_PARAMFILE() 
;
; PURPOSE:
;   Read a SFHGRID parameter file.
;
; INPUTS: 
;   sfhgrid_paramfile - parameter file to read
;
; OPTIONAL INPUTS: 
;   thissfhgrid - subscript the file to this SFH grid (default is to
;     read them all)
;
; OUTPUTS: 
;   allparams - data structure with all the parameters 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Nov 12, UCSD
;   jm13jan13siena - SFHGRID_PARAMFILE now required input
;
; Copyright (C) 2010, 2013, John Moustakas
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

function read_sfhgrid_paramfile, sfhgrid_paramfile, thissfhgrid=thissfhgrid

    if (n_elements(sfhgrid_paramfile) eq 0) then message, $
      'SFHGRID_PARAMFILE input required!'
    if (file_test(sfhgrid_paramfile,/regular) eq 0) then $
      message, 'SFHGRID parameter file '+sfhgrid_paramfile+' not found'

    allparams = yanny_readone(sfhgrid_paramfile)

; build the parameter arrays for the specified SFHGRID
    if (n_elements(thissfhgrid) ne 0) then begin
       match = where(allparams.sfhgrid eq thissfhgrid)
       if (match[0] eq -1) then message, 'Please update '+sfhgrid_paramfile
       params = allparams[match]
       return, params
    endif

return, allparams
end
    
