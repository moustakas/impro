;+
; NAME:
;   READ_ISEDFIT_PARAMFILE()
;
; PURPOSE:
;   Read an iSEDfit-style parameter file.
;
; INPUTS: 
;   paramfile - parameter file name
;
; OPTIONAL INPUTS: 
;   sfhgrid - if sfhgrid in the parameter file is an array then
;     subscript to SFHGRID
;
; OUTPUTS: 
;   params - data structure containing the specified parameters 
;
; COMMENTS:
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

function read_isedfit_paramfile, paramfile, sfhgrid=sfhgrid

    if (file_test(paramfile,/regular) eq 0) then $
      message, 'PARAMFILE '+paramfile+' not found'

    lines1 = djs_readlines(paramfile)
    keep = where((strcompress(lines1,/remove) ne '') and $
      (strmatch(lines1,'#*') eq 0),nkeep)
    if (nkeep eq 0) then message, 'Problem parsing parameter file '+$
      paramfile
    lines = lines1[keep]

    for ii = 0, n_elements(lines)-1 do begin
       words = strsplit(lines[ii],' ',/extract)
       words = words[where(strcompress(words,/remove) ne '')]
       name = strtrim(words[0],2)
       value1 = strtrim(words[1],2)
       case name of
          'sfhgrid': begin
             if (n_elements(sfhgrid) eq 0) then begin
                value = fix(strsplit(value1,',',/extract))
                if (n_elements(value) eq 1) then value = value[0]
             endif else value = sfhgrid ; overwrite the parameter file
          end
          'redcurve': begin
             value = strsplit(value1,',',/extract)
             if (n_elements(value) eq 1) then value = value[0]
          end
          'filterlist': value = strsplit(value1,',',/extract)
          'igm': value = fix(value1)
          'maxold': value = fix(value1)
          'redshift': begin
             zsplit = im_double(strsplit(value1,',',/extract))
             minz = zsplit[0]
             maxz = zsplit[1]
             nz = zsplit[2]
             if (nz le 1) then message, 'Update your parameter file to use NZ instead of dz!!'
             value = range(minz,maxz,nz)
;            value = im_array(minz,maxz,dz,/double)
          end
          else: value = value1
       endcase
       if (ii eq 0) then params = create_struct(name,value) else $
         params = create_struct(params,name,value)
    endfor

; add any missing keywords here    
    if (tag_exist(params,'maxold') eq 0) then params = $
      struct_addtags(params,{maxold: 0})
    
;; if SFHGRID is a vector, then change PARAMS into an array of
;; structures
;    nsfhgrid = n_elements(params.sfhgrid)
;    if (nsfhgrid gt 1) then begin
;       tags = tag_names(params)
;       for ii = 0, n_elements(tags)-1 do begin
;          if (ii eq 0) then $
;            newparams = create_struct(tags[ii],(params.(ii))[0]) else $
;              newparams = create_struct(newparams,tags[ii],(params.(ii))[0])
;       endfor
;    endif

return, params
end
    
