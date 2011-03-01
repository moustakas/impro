pro im_mwrfits, outstruct, outfile, hdr, silent=silent, $
  append=append, gzip=gzip, nogzip=nogzip, clobber=clobber
; jm09mar18nyu - simple wrapper on MWRFITS that does what I want
; jm09nov10ucsd - added /CLOBBER keyword

    if keyword_set(append) then begin
       if (file_test(outfile) eq 0) then begin
          splog, 'Output file '+outfile+' does not exist!'
          return
       endif
       if (keyword_set(silent) eq 0) then splog, 'Appending to '+outfile
       mwrfits, outstruct, outfile, hdr, /silent
       if keyword_set(gzip) then spawn, 'gzip -f '+outfile, /sh
    endif else begin
       if im_file_test(repstr(outfile,'.gz','')+'.gz',clobber=clobber) then return
       if (keyword_set(silent) eq 0) then splog, 'Writing '+outfile
       mwrfits, outstruct, outfile, hdr, /create
       if (keyword_set(nogzip) eq 0) then spawn, 'gzip -f '+outfile, /sh
    endelse
    
return
end
    
