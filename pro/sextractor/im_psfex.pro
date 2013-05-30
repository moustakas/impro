;+
; NAME:
;   IM_PSFEX
;
; PURPOSE:
;   Run PSFEx on a suitably prepared SExtractor catalog.
;
; INPUTS: 
;   sexcatalog - name of input SE catalog
;
; OPTIONAL INPUTS: 
;   config - PSFEx configuration structure (see INIT_PSFEX_CONFIG); if
;     none is provided then a default configuration file is used
;   configfile - optional output configuration file name for
;     WRITE_DUMMY_CONFIG 
;
; KEYWORD PARAMETERS: 
;   silent - suppress messages to STDOUT
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
;   J. Moustakas, 2013 May 24, Siena
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

pro im_psfex, sexcatalog, config, configfile=configfile, silent=silent

    ncat = n_elements(sexcatalog)
    nconfig = n_elements(config)

    if (ncat eq 0L) then begin
       doc_library, 'im_psfex'
       return
    endif

; make a default CONFIG structure    
    if (nconfig eq 0L) then begin
       config = init_psfex_config(ncat)
       nconfig = ncat
    endif
    
    if (ncat ne nconfig) then begin
       splog, 'Dimensions of SEXCATALOG and CONFIG do not match'
       return
    endif

    for ii = 0L, ncat-1L do begin
       if (not keyword_set(silent)) then splog, $
         'Running PSFEx '+sexcatalog[ii]
       write_dummy_config, config[ii], configfile=configfile
       spawn, 'psfex '+sexcatalog[ii]+' -c '+configfile
    endfor
    
return
end
