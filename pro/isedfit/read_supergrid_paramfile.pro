;+
; NAME:
;   READ_SUPERGRID_PARAMFILE() 
;
; PURPOSE:
;   Read a SUPERGRID parameter file.
;
; INPUTS: 
;   supergrid_paramfile - parameter file to read
;
; OPTIONAL INPUTS: 
;   supergrid - subscript the file to this SUPERGRID (default is to read
;     them all)
;
; OUTPUTS: 
;   allparams - data structure with all the parameters 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Nov 12, UCSD
;   jm13jan13siena - SUPERGRID_PARAMFILE now required input
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

function read_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid

    if (n_elements(supergrid_paramfile) eq 0) then message, $
      'SUPERGRID_PARAMFILE input required!'
    if (file_test(supergrid_paramfile,/regular) eq 0) then $
      message, 'SUPERGRID parameter file '+supergrid_paramfile+' not found'

    allparams = yanny_readone(supergrid_paramfile)

; build the parameter arrays for the specified SUPERGRID
    if (n_elements(supergrid) ne 0) then begin
       match = where(allparams.supergrid eq supergrid)
       if (match[0] eq -1) then message, 'No matching SUPERGRID '+$
         strtrim(supergrid,2)+' in parameter file '+supergrid_paramfile
       params = allparams[match]
       return, params
    endif

return, allparams
end
    
