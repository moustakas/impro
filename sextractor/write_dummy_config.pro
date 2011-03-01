;+
; NAME:
;   WRITE_DUMMY_CONFIG
;
; PURPOSE:
;   Write a dummy (temporary) configuration file (for use with
;   IM_SEX, IM_SCAMP, and IM_SWARP).
;
; INPUTS: 
;   config - configuration structure
;
; OPTIONAL INPUTS: 
;   configfile - output file name (default is to create a dummy
;   file in '/tmp/')
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
;   J. Moustakas, 2008 Aug 06, NYU - written
;
; Copyright (C) 2008, John Moustakas
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

pro write_dummy_config, config, configfile=configfile

    nconfig = n_elements(config)
    if (nconfig eq 0L) then begin
       doc_library, 'write_dummy_config'
       return
    endif

    if (nconfig gt 1L) then begin
       if (n_elements(configfile) eq 0L) then begin
          configfile = '/tmp/dummy.'+strtrim(lindgen(nconfig),2)+'.config'
       endif else begin
          if (nconfig ne n_elements(configfile)) then begin
             splog, 'Dimensions of CONFIG and CONFIGFILE must match'
             return
          endif
       endelse
       for ic = 0L, nconfig-1L do $
         write_dummy_config, config[ic], configfile[ic]
       return
    endif
    
    if (n_elements(configfile) eq 0L) then $
      configfile = '/tmp/dummy.config'

    tags = tag_names(config[0])
    
    openw, lun, configfile, /get_lun
    for ii = 0L, nconfig-1L do $
      for jj = 0L, n_tags(config)-1L do $
        printf, lun, tags[jj]+' '+strtrim(config[ii].(jj),2)
    free_lun, lun
        
return
end
    
