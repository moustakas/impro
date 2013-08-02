;+
; NAME:
;   READ_ISEDFIT()
;
; PURPOSE:
;   Simple function to read the output of iSEDfit.
;
; INPUTS:
;   paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   isedfit_dir - I/O path
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   isedfit - ISEDFIT result structure
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 17, Siena
;
; Copyright (C) 2013, John Moustakas
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

function read_isedfit, isedfit_paramfile, params=params, index=index, $
  supergrid_paramfile=supergrid_paramfile, thissupergrid=thissupergrid, $
  isedfit_dir=isedfit_dir

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'read_isedfit'
       return, -1
    endif

    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile)

; read the SUPERGRID parameter file
    if n_elements(supergrid_paramfile) eq 0 then begin
       splog, 'SUPERGRID parameter file required'
       return, -1
    endif
    
    super = read_supergrid_paramfile(supergrid_paramfile,supergrid=thissupergrid)
    if n_elements(thissupergrid) eq 0 then thissupergrid = super.supergrid

    fp = isedfit_filepaths(params,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=thissupergrid,isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
      outprefix=outprefix)
    isedfile = fp.isedfit_dir+fp.isedfit_outfile+'.gz'
    if (file_test(isedfile,/regular) eq 0) then begin
       splog, 'iSEDfit file '+isedfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+isedfile
    isedfit = mrdfits(isedfile,1,/silent,rows=index)
    
return, isedfit
end
