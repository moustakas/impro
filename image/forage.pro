
; returns path name of the directory your data files are in. 
; called by 'forage.pro'

FUNCTION datapath

; the object frames

  pathname = '/deep1/ioannis/chile/ctio/n1obj/'
  return, pathname

END

; extract header information from FITS header h, and put it into 
; a dummy structure called 'a'
; called by 'forage.pro' for every image

FUNCTION head_extract, path, fname, template

  fxread, path+fname, dum, h, 1, 2, 1, 2
  a = template
  a.fname   = fname
  a.object  = fxpar(h, 'object')
  a.filter  = fxpar(h, 'filter1')
  a.type    = fxpar(h, 'imagetyp')
  a.exptime = fxpar(h, 'exptime')
  a.utstart = fxpar(h, 'utshut')
  a.ra      = fxpar(h, 'ra')
  a.dec     = fxpar(h, 'dec')

  return, a
END 


; prints the information found by 'forage.pro' to the screen
; called by 'forage.pro'

PRO printinfo, gen, frame, filter=filter, abbrev=abbrev, contains=contains

  print
  print, 'Telescope:     ', gen.telescope
  print, 'Detector:      ', gen.detector, $
    string(gen.xsize, format='(I7)'), ' x', string(gen.xsize, format='(I5)')
  print, 'Observered by: ', gen.observer
  print, replicate('-', 39)
  print, '     Filename           Object    Filter    Type     Exptime    UTStart'
  print, replicate('-', 39)

  nframes = n_elements(frame)
  FOR i=0, nframes-1 DO BEGIN 
     fi = frame[i]
     IF keyword_set(abbrev) THEN BEGIN 
; to abbreviate the list, only list those frames in filter 'i'
        IF fi.filter EQ 'i' THEN BEGIN 
           print, fi.fname, fi.object, fi.filter, fi.type, fi.exptime, $
            fi.utstart, format='(A15,A20,A3,A15,F7.1,A14)'
        ENDIF 
     ENDIF ELSE BEGIN 
; only print out those frames in a specified filter
;        IF keyword_set(filter) THEN BEGIN 
;           dstrng = string(fi.fname, fi.object, fi.filter, fi.type, $
;                           fi.exptime, fi.utstart, $
;                           format='(A15,A20,A3,A15,F7.1,A14)')
;           IF (strpos(fi.filter,'filter') GT -1 THEN BEGIN 
;              print, dstrng
;           ENDIF 
 ;       ENDIF ELSE BEGIN 
        IF keyword_set(contains) THEN BEGIN 
; search for a particular string in the structure
           dstrng = string(fi.fname, fi.object, fi.filter, fi.type, $
                           fi.exptime, fi.utstart, $
                           format='(A15,A20,A3,A15,F7.1,A14)')
           IF  (strpos(strupcase(dstrng),strupcase(contains))) GT -1 THEN BEGIN
              print, dstrng
           ENDIF 
        ENDIF ELSE BEGIN
           print, fi.fname, fi.object, fi.filter, fi.type, $
            fi.exptime, fi.utstart, format='(A15,A20,A3,A15,F7.1,A14)'
        ENDELSE 
     ENDELSE 
  ENDFOR 
  
  print
  print, nframes, ' frames'
  return
END




; create general information and image-specific structures

PRO forage, gen, frame, silent=silent
   
   gen = {gen_info, $
          telescope: '', $
          detector:  '', $
          observer:  '', $
          nampsyx:   fltarr(2, 2), $
          xsize:     0, $
          ysize:     0}  
   
   frame_template = {frame_info, $
                     fname:   '', $
                     object:  '', $
                     filter:  '', $
                     type:    '', $
                     exptime: 0.0, $
                     utstart: '', $
                     ra:    '', $
                     dec:   ''}
   
   flist = findfile(datapath())

; keep only fits files

  fgood = where(strpos(flist, '.fits') GE 0)
  flist = flist[fgood]
  fcount = n_elements(flist)

  IF fcount GE 0 THEN BEGIN 

; read little 2x2 subimage in the corner to save time

      fxread, datapath()+flist[0], dum, h, 1, 2, 1, 2
  ENDIF ELSE BEGIN 
      print, 'We want FITS files!!'
      return
  ENDELSE 

; create an array of structures for the total number of images

  frame = replicate(frame_template, fcount)

  gen.telescope = fxpar(h, 'telescop') 
  gen.detector  = fxpar(h, 'detector')
  gen.observer  = fxpar(h, 'observer')
  gen.nampsyx   = fxpar(h, 'nampsyx')
  gen.xsize     = fxpar(h, 'naxis1')
  gen.ysize     = fxpar(h, 'naxis2')

  FOR i=0, fcount-1 DO frame[i] = head_extract(datapath(), flist[i], frame_template)

; call 'printinfo' and print everything to screen

  IF keyword_set(silent) EQ 0 THEN printinfo, gen, frame

  return
END
