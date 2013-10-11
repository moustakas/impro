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
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     read this SFHgrid (may be a vector)
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

; build the redshift array; WRITE_ISEDFIT_PARAMFILE makes sure that
; NZZ is the same for each SFHGRID
    params = struct_addtags(params,replicate({redshift: $
      fltarr(params[0].nzz)},n_elements(params)))
    for ii = 0, n_elements(params)-1 do begin
       if params[ii].user_redshift then params[ii].redshift = params[ii].use_redshift else begin
          if params[ii].nzz gt 1 then begin
             params[ii].redshift = range(params[ii].zminmax[0],params[ii].zminmax[1],$
               params[ii].nzz,log=params[ii].zlog)
          endif else begin
             params[ii].redshift = params[ii].zminmax[0] ; assume ZMIN=ZMAX
          endelse
       endelse
    endfor
    
; add NMAXBURST and NMODELCHUNK (see ISEDFIT_MONTEGRIDS)
    params = struct_addtags(params,replicate({nmaxburst: 0L, $
      nmodelchunk: 0L},n_elements(params)))
    for ii = 0, n_elements(params)-1 do begin
       if (params[ii].pburst gt 0D) then params[ii].nmaxburst = $
         ceil((params[ii].tburst[1]-params[ii].tburst[0])/params[ii].interval_pburst)
       params[ii].nmodelchunk = ceil(params[ii].nmodel/float(params[ii].modelchunksize))
    endfor
    
; choose one particular SFHGRID
    if (n_elements(thissfhgrid) ne 0) then begin
       match = where(params.sfhgrid eq thissfhgrid)
       if (match[0] eq -1) then message, 'No matching SFHgrid '+$
         strtrim(thissfhgrid,2)+'; please update '+isedfit_paramfile
       params = params[match]
       return, params
    endif

return, params
end
    
