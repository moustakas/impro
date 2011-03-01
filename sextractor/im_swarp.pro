;+
; NAME:
;       IM_SWARP
;
; PURPOSE:
;       Run SWARP given a list of images.
;
; INPUTS: 
;       imagelist - list of FITS file names
;       config - SWARP configuration structure (see
;         INIT_SWARP_CONFIG); must be a scalar structure 
;
; OPTIONAL INPUTS: 
;       configfile - optional output configuration file name for
;         WRITE_DUMMY_CONFIG 
;
; KEYWORD PARAMETERS: 
;       silent - suppress messages to STDOUT
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
;       J. Moustakas, 2008 Aug 06, NYU - written
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

pro im_swarp, imagelist, config, configfile=configfile, silent=silent

    nimage = n_elements(imagelist)
    nconfig = n_elements(config)

    if (nimage eq 0L) or (nconfig eq 0L) then begin
       doc_library, 'im_swarp'
       return
    endif

    if (nconfig ne 1L) then begin
       splog, 'CONFIG must be a scalar structure'
       return
    endif

    null = where(strtrim(imagelist,2) eq '',nnull)
    if (nnull ne 0) then begin
       splog, 'At least some elements of IMAGELIST are blank'
       return
    endif
    
;   if (not keyword_set(silent)) then splog, 'Spawning SWARP'
    write_dummy_config, config, configfile=configfile
    spawn, 'swarp '+strjoin(imagelist,',')+' -c '+configfile
    
return
end
