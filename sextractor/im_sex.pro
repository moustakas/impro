;+
; NAME:
;   IM_SEX
;
; PURPOSE:
;   Run SExtractor on one or more images (also optionally in
;   double-image mode). 
;
; INPUTS: 
;   imagelist - list of FITS file names
;
; OPTIONAL INPUTS: 
;   config - SE configuration structure (see INIT_SEX_CONFIG); if none
;     is provided then a default configuration file is used
;   configfile - optional output configuration file name for
;     WRITE_DUMMY_CONFIG 
;   detect_imagelist - image(s) to use for detection (double-image 
;     mode); can be a vector or a scalar
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
;   J. Moustakas, 2001 Aug 30, U of A - written
;   jm08aug06nyu - major overhaul
;
; Copyright (C) 2001, 2008, John Moustakas
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

pro im_sex, imagelist, config, detect_imagelist=detect_imagelist1, $
  configfile=configfile, silent=silent

    nimage = n_elements(imagelist)
    nconfig = n_elements(config)

    if (nimage eq 0L) then begin
       doc_library, 'im_sex'
       return
    endif

; make a default CONFIG structure    
    if (nconfig eq 0L) then begin
       config = init_sex_config(nimage)
       nconfig = nimage
    endif
    
    if (nimage ne nconfig) then begin
       splog, 'Dimensions of IMAGELIST and CONFIG do not match'
       return
    endif

    if (n_elements(detect_imagelist1) ne 0L) then begin ; double-image mode
       if (n_elements(detect_imagelist1) eq 1L) then $
         detect_imagelist = replicate(detect_imagelist1,nimage) else begin
          if (n_elements(detect_imagelist1) ne nimage) then begin
             splog, 'Dimensions of IMAGELIST and DETECT_IMAGELIST do not match'
             return
          endif else detect_imagelist = detect_imagelist1
       endelse
       double = 1
    endif
    
    for ii = 0L, nimage-1L do begin
       if (not keyword_set(silent)) then splog, $
         'SExtracting '+imagelist[ii]
       write_dummy_config, config[ii], configfile=configfile
       if keyword_set(double) then begin
          spawn, 'sex '+detect_imagelist[ii]+','+imagelist[ii]+' -c '+configfile
       endif else begin
          spawn, 'sex '+imagelist[ii]+' -c '+configfile
       endelse
    endfor
    
return
end
