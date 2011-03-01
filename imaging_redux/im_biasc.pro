

;------------------------------------------------------------------
;+
; NAME:
;    rjc_biasc
;
; PURPOSE:
;    Bias subtract given image
;
;  CALLING SEQUENCE:
;      rjc_biasc, image, biasframe=biasframe, hdr=hdr
;
;  INPUTS:
; 
;     image : [nx, ny] image to be overscan corrected
;     biasframe : [nx, ny] bias image
;
;  OPTIONAL INPUTS:
;     
;     
;  -
;-------------------------------------------------------------

pro im_biasc, image, biasframe=biasframe, biasfile=biasfile, hdr=hdr, silent=silent
   
   nx = n_elements(image[*,0])
   ny = n_elements(image[0,*])
   
   bnx = n_elements(biasframe[*,0])
   bny = n_elements(biasframe[0,*])
   
   IF bnx NE nx OR bny NE ny THEN BEGIN
      splog, 'The dimensions of the biasframe and input image do not match'
      stop
   ENDIF

   if (n_elements(biasfile) eq 0L) then biasfile = 'Unknown'

   if (not keyword_set(silent)) then splog, 'Subtracting the bias frame, '+biasfile
   
   ;;Perform the subtraction
   image = image - biasframe
   
   if keyword_set(hdr) then begin
      sxaddpar, hdr, 'BIASCOR', biasfile, ' bias image', before='HISTORY'
      sxaddhist, "'Bias frame subtracted "+hogg_iso_date()+"'", hdr
   endif
   
return
end
