;+
; NAME:
;   IM_FILE_TEST()
; PURPOSE:
;   Check if a file exists before overwriting, unless /CLOBBER. 
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Sep 12, UCSD
;-

function im_file_test, file, clobber=clobber
    if file_test(file) and (keyword_set(clobber) eq 0) then begin
       splog, 'File '+file+' exists; use /CLOBBER'
       return, 1
    endif
return, 0
end
