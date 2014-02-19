;+
;  im_update_header() -- updates the header with units, and comments.
;  Input:
;     filename -- path to fits file
;     units -- array of "tag,unit,comment" strings -- see example
;     nheader -- starting array of tags
;     /write -- update the filename
;  Returns:
;     header array.
;  History:
;     Mendez 2014 from primus_website.
;  
;-


FUNCTION im_update_header, filename, units, nheader, write=write
  ;; Fix up the units, this is a quick an dirty hack...
  ;; units = ["tag,unit,comment", ...]
  ;; nheader = ["X=item", "X=item", ...]

  cat =  mrdfits(filename, 1, header)
  nl =  string(10B) ; new line
  
  forprint, header

  ;; do a quick loop to add units to data header
  FOR j=0, n_elements(header)-1 DO BEGIN
    modified = 0
    FOR i=0, n_elements(units)-1 DO BEGIN
      tmp = strsplit(units[i], ',', /extract,  /preserve_null)
      tag = strupcase(strtrim(tmp[0], 2))
      unit = tmp[1]
      comment = tmp[2]

      IF strmatch(header[j], 'TTYPE*') THEN BEGIN
        ;; get the column name
        tmp = strsplit(header[j], "'", /extract,  /preserve_null)
        col = strupcase(strtrim(tmp[1], 2))
        
        IF strmatch(tag, col) THEN BEGIN
          n = strmid(header[j], 5, 2)
          nheader = [nheader, $
                     ; add in the column name with coment
                     (strsplit(header[j], '/', /extract))[0]+'/ '+comment,  $
                     ; add in the units -- python picks this out
                     'TUNIT'+n+" = '"+string(unit, '(A-8)')+"'           /"]
          modified = 1
        ENDIF
      ENDIF
    ENDFOR 
    IF modified EQ 0 AND strmatch(header[j], 'TTYPE*') THEN BEGIN
      nheader = [nheader, header[j]]
    ENDIF
  ENDFOR

  header = [nheader, "END"]
  
  ;; Need to strip type creation so that I can populate the comments
  IF keyword_set(write) THEN BEGIN
    mwrfits, cat, repstr(filename, '.gz', ''), header, /create,  /no_type
    IF strmatch(filename, '*.gz') THEN $
       spawn, 'gzip -f '+repstr(filename, '.gz', ''), /sh
    splog, 'Updated Units: ', filename
  ENDIF
  
  return, header
END

PRO update_header_example
  filename = 'test.fits'
  mwrfits, {ra:0, dec:0, z:0,  zwarn:0}, filename,  /create
  
  
  ;; array of unit, and commets for the columns of the fits file.
  units = ['RA,degree,Right Ascension J2000', $
           'DEC,degree,Declination J2000', $
           'Z,,Best-fit PRIMUS redshift']
  
  ;; starting set of tags
  nheader =  ["SURVEY    = 'PRIsm MUlti-Object Survey (PRIMUS)'", $
              "DATASET   = 'Redshift Catalog'", $
              "VERSION   = '2.0'", $
              "FILENAME  = '"+filename+"'", $
              "DATE      = '"+systime()+"'", $
              "PAPER1    = '2011ApJ,741,8C'", $
              "PAPER1URL = 'http://adsabs.harvard.edu/abs/2011ApJ...741....8C'",$
              "PAPER2    = '2013ApJ,767,118C'", $
              "PAPER2URL = 'http://adsabs.harvard.edu/abs/2013ApJ...767..118C'", $
              "COMPILEDBY= 'Alexander Mendez (ajmendez@physics.ucsd.edu)'", $
              ;; "WRITTENBY = '"+(reverse(scope_traceback()))[0]+"'", $
              "URL       = 'http://primus.ucsd.edu/'   / For more info about the PRIMUS dataset"]

  header =  im_update_header(filename, units, nheader,  /write)
  
  x = mrdfits(filename, 1, hdr)
  forprint, hdr
END
