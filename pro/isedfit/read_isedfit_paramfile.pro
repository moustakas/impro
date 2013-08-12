;+
; NAME:
;   READ_ISEDFIT_PARAMFILE() 
;
; PURPOSE:
;   Read the iSEDFIT parameter file.
;
; INPUTS: 
;   isedfit_paramfile - parameter file to read
;
; OPTIONAL INPUTS: 
;   thissfhgrid - subscript the file to this SFH grid (default is to
;     read them all)
;
; OUTPUTS: 
;   params - data structure with all the parameters 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Nov 12, UCSD
;   jm13jan13siena - ISEDFIT_PARAMFILE now required input
;   jm13aug05siena - updated to latest data model
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

function read_isedfit_paramfile, isedfit_paramfile, thissfhgrid=thissfhgrid

    if n_elements(isedfit_paramfile) eq 0 then message, $
      'ISEDFIT_PARAMFILE input required!'
    if file_test(isedfit_paramfile) eq 0 then $
      message, 'SFHGRID parameter file '+isedfit_paramfile+' not found'

    params = yanny_readone(isedfit_paramfile)

; add NMAXBURST (see ISEDFIT_BUILD_MONTEGRIDS)
    params = struct_addtags(params,replicate({nmaxburst: 0L},n_elements(params)))
    for ii = 0, n_elements(params)-1 do if (params[ii].pburst gt 0D) then $
      params[ii].nmaxburst = ceil((params[ii].tburst[1]-params[ii].tburst[0])/$
      params[ii].interval_pburst)
    
; build the parameter arrays for the specified SFHGRID
    if (n_elements(thissfhgrid) ne 0) then begin
       match = where(params.sfhgrid eq thissfhgrid)
       if (match[0] eq -1) then message, 'No matching SFHgrid '+$
         strtrim(thissfhgrid,2)+'; please update '+isedfit_paramfile
       params = params[match]
       return, params
    endif

return, params
end
    
