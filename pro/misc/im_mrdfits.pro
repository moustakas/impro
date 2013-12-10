;+
; NAME:
;   im_mrdfits
; PURPOSE:
;   wrapper around mrdfits to grab a possible fits file
; CALLING SEQUENCE:
;   str= gz_mrdfits(filename, [ext], [hdr], [mrdfits extra])
; COMMENTS:
;   pass in a exact file, or just the filename excluding the extension
; REVISION HISTORY:
;   2013.12.04 -- Mendez
;-
;------------------------------------------------------------------------------
FUNCTION im_mrdfits, filename, ext, hdr, $
                     status=status, silent=silent, $
                     _extra=extra
  IF n_elements(ext) EQ 0 THEN ext = 1
  
  IF file_test(filename) EQ 0 THEN BEGIN
    ;; could not find the file, so search for similar files
    ;; that have an additional ending -- compressed things
    files = file_search(filename+'*', /fold_case)
    IF n_elements(files) NE 1 THEN BEGIN
      print, 'Found the following possible files: '
      forprint, files, format='("     ",(A))'
      message, 'Could not uniqly determine file to open, please specify file better.'
    ENDIF ELSE BEGIN
      filename = files[0]
    ENDELSE
  ENDIF
  IF NOT keyword_set(silent) THEN splog, 'Loading: '+filename
  return, mrdfits(filename, ext, hdr, status=status, silent=silent, _extra=extra)
END

  
  
