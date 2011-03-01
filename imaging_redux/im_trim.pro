;------------------------------------------------------------------
;+
; NAME:
;    rjc_trim
;
; PURPOSE:
;    trim input image
;    
;
;
;  CALLING SEQUENCE:
;
;      rjc_trim, image, trimsec, hdr=hdr
;
;  INPUTS:
; 
;     image : [nx, ny] image to be overscan corrected
;     trimsec : section to use for trim [x0, x1, y0, y1]
;
;  OPTIONAL INPUTS:
;     
;  -
;-------------------------------------------------------------
pro im_trim, image, trimsec=trimsec, hdr=hdr, silent=silent
   
   IF NOT keyword_set(trimsec) THEN BEGIN
      IF NOT keyword_set(hdr) THEN BEGIN
         splog, 'You must specify either a hdr or a trimsec'
         stop
      ENDIF ELSE BEGIN
         bs = sxpar(hdr, 'TRIMSEC')
         trimsec = strsplit(bs, '[]:,', /extract)
         ;;Now correct for the fact that IRAF uses 1 indexing and
         ;;IDL uses zero
         trimsec = trimsec - 1
      ENDELSE
   ENDIF

   IF n_elements(image) EQ 0 THEN BEGIN
      splog, 'Please input an image with more than 0 pixels'
      stop
   endIF
   
   if (not keyword_set(silent)) then splog, 'Trimming'

   nx = n_elements(image[*,0])
   ny = n_elements(image[0,*])
   
   outputimage = image[trimsec[0]:trimsec[1],trimsec[2]:trimsec[3]]
   
   image = outputimage
   
; if the hdr was given to the program, update it
   
   IF keyword_set(hdr) THEN begin
      sxaddpar, hdr, 'TRIM', '['+string(trimsec[0]+1,format='(I0)')+':'+string(trimsec[1]+1,format='(I0)')+','+$
        string(trimsec[2]+1,format='(I0)')+':'+string(trimsec[3]+1,format='(I0)')+']', ' trim region', before='HISTORY'
      sxaddhist, "'Image trimmed "+hogg_iso_date()+"'", hdr
   ENDIF 
   
return
end
