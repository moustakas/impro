;------------------------------------------------------------------
;+
; NAME:
;    rjc_darkc
;
; PURPOSE:
;    Dark subtract given image
;
;  CALLING SEQUENCE:
;      rjc_darkc, image, darkframe=darkframe, hdr=hdr
;
;  INPUTS:
; 
;     image : [nx, ny] image to be overscan corrected
;     darkframe : [nx, ny] dark image [counts/s]
;
;  OPTIONAL INPUTS:
;     
;     
;  -
;-------------------------------------------------------------
pro im_darkc, image, exptime, darkframe=darkframe, darkfile=darkfile, $
  hdr=hdr, silent=silent
   
   nx = n_elements(image[*,0])
   ny = n_elements(image[0,*])
   
   bnx = n_elements(darkframe[*,0])
   bny = n_elements(darkframe[0,*])
   
   IF bnx NE nx OR bny NE ny THEN BEGIN
      splog, 'The dimensions of the darkframe and input image do not match'
      stop
   ENDIF
   
   if (n_elements(darkfile) eq 0L) then darkfile = 'Unknown'

   if (not keyword_set(silent)) then splog, 'Subtracting the dark frame, '+darkfile+'.'

   image = image - exptime*darkframe ; darkframe must be in counts/s

   ;;Update the hdr

   if keyword_set(hdr) then begin
      sxaddpar, hdr, 'DARKCOR', darkfile, ' dark image', before='HISTORY'
      sxaddhist, "'Dark frame subtracted "+hogg_iso_date()+"'", hdr
   endif
   
return
end
