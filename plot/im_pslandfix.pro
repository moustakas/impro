pro im_pslandfix, filename, xw=xw, yw=yw
;+
; NAME:
;   PSLANDFIX
;
; PURPOSE:
;   Fix the upside-down landscape that IDL generates.
;
; CALLING SEQEUNCE:
;   pro pslandfix,[filename]
;
; INPUT:
;   FILENAME   Name the idl PostScript file to fix.  Default is 'idl.ps'
;
; OUTPUT:
;   Fixed PosctScript file
; NOTES:
;   none
; HISTORY:
;   10-AUG-95 Version 1 written    E. Deutsch
;-

  psfile = djs_readlines(filename)
  these = where(strmatch(psfile,'*$IDL_DICT*') and $
    strmatch(psfile,'*rotate*'),nthese)
  for ii = 0L, nthese-1L do begin
     x1 = strpos(psfile[these[ii]],'scale')
     x2 = strpos(psfile[these[ii]],'rotate')
     doit = psfile[these[ii]]
     strput, doit,' 90', x1+6
     doit = repstr(doit,'0 792 ','560 40 ')
     psfile[these[ii]] = doit
  endfor

  test = repstr(filename,'.ps','.test.ps')
  spawn, '/bin/cp '+filename+' '+test, /sh
  openw, lun, test, /get_lun
  for jj = 0L, n_elements(psfile)-1L do $
    printf, lun, psfile[jj]
  free_lun, lun
  
return
end
