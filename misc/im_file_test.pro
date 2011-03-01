function im_file_test, file, clobber=clobber
; jm10sep12ucsd - check if a file exists before overwriting, unless
; /CLOBBER 
    if file_test(file) and (keyword_set(clobber) eq 0) then begin
       splog, 'File '+file+' exists; use /CLOBBER'
       return, 1
    endif
return, 0
end
