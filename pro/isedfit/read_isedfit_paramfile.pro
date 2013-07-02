;+
; NAME:
;   READ_ISEDFIT_PARAMFILE()
;
; PURPOSE:
;   Read an iSEDfit-style parameter file.
;
; INPUTS: 
;   isedfit_paramfile - parameter file name
;
; OPTIONAL INPUTS: 
;   use_redshift - use this redshift array instead of constructing the
;     redshift array from the parameters given in the
;     ISEDFIT_PARAMFILE parameter file; useful for when you have a
;     relatively sample of objects with well-determined redshifts
;     spanning a wide redshift range [NZZ]
;
; OUTPUTS: 
;   params - data structure containing the specified parameters 
;
; COMMENTS:
;   Needs better error checking to ensure the parameter file is in the
;   right format.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 11, NYU - written
;   jm09may16nyu - changed the parameter file format
;   jm10feb11ucsd - documented
;
; Copyright (C) 2009-2010, John Moustakas
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

function read_isedfit_paramfile, isedfit_paramfile, use_redshift=use_redshift

    if (file_test(isedfit_paramfile,/regular) eq 0) then $
      message, 'PARAMFILE '+isedfit_paramfile+' not found'

    lines1 = djs_readlines(isedfit_paramfile)
    keep = where((strcompress(lines1,/remove) ne '') and $
      (strmatch(lines1,'#*') eq 0),nkeep)
    if (nkeep eq 0) then message, 'Problem parsing parameter file '+$
      isedfit_paramfile
    lines = lines1[keep]

    for ii = 0, n_elements(lines)-1 do begin
       words = strsplit(lines[ii],' ',/extract)
       words = words[where(strcompress(words,/remove) ne '')]
       name = strtrim(words[0],2)
       if (n_elements(words) eq 1) then value1 = '' else $ ; blank, probably REDCURVE
         value1 = strtrim(words[1],2)
       case name of
          'filterlist': value = strsplit(value1,',',/extract)
          'igm': value = fix(value1)
          'maxold': value = fix(value1)
          'minz': value = double(value1)
          'maxz': value = double(value1)
          'nzz': value = long(value1)
          'zlog': value = fix(value1)
          'h100': value = double(value1)
          'omega0': value = float(value1)
          'omegal': value = float(value1)
          else: value = value1
       endcase
       if (ii eq 0) then params = create_struct(name,value) else $
         params = create_struct(params,name,value)
    endfor

; add any missing keywords here    
    if (tag_exist(params,'maxold') eq 0) then params = $
      struct_addtags(params,{maxold: 0})
    
; build the redshift array, with the option of overwriting the default
; parameters 
    nzz = n_elements(use_redshift)
    if nzz ne 0 then begin
       params.nzz = nzz
       params.minz = min(use_redshift)
       params.maxz = max(use_redshift)
       params = struct_addtags(params,{redshift: use_redshift})
    endif else begin
       if params.nzz eq 1 then params = struct_addtags(params,{redshift: 0D}) else $
         params = struct_addtags(params,{redshift: dblarr(params.nzz)})
       params.redshift = range(params.minz,params.maxz,params.nzz,log=params.zlog eq 1)
    endelse

return, params
end
