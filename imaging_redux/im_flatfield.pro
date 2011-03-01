;+
; NAME:
;    im_flatfield
;
; PURPOSE:
;    Flat-field a given image.
;
;  CALLING SEQUENCE:
;      im_flatfield, 
;
;  INPUTS:
;     image : [nx, ny] image to be overscan corrected
;     biasframe : [nx, ny] bias image
;
;  OPTIONAL INPUTS:
;     
;  -

pro im_flatfield, image, varmap=varmap, flatfield=flatfield, $
  flatfile=flatfile, hdr=hdr, silent=silent
   
   nx = n_elements(image[*,0])
   ny = n_elements(image[0,*])
   
   fnx = n_elements(flatfield[*,0])
   fny = n_elements(flatfield[0,*])
   
   if (fnx ne nx) or (fny ne ny) then begin
      splog, 'Dimensions of the flatfield and input image do not match'
      return
   endif

   if arg_present(varmap) then begin
      if (n_elements(varmap) ne n_elements(image)) then begin
         splog, 'Dimensions of IMAGE and VARMAP do not match'
         return
      endif
   endif
   
   if (n_elements(flatfile) eq 0L) then flatfile = 'unknown'
   if (not keyword_set(silent)) then splog, 'Dividing by the flat-field, '+flatfile

; divide by the flat-field and update the header; protect against
; divide-by-zero 

   image = image / (flatfield + (flatfield eq 0.0))
   if arg_present(varmap) then varmap = varmap / (flatfield + (flatfield eq 0.0))^2.0

   if keyword_set(hdr) then begin
      sxaddpar, hdr, 'FLATCOR', flatfile, ' flat-field', before='HISTORY'
      sxaddhist, "'Flat-fielded "+hogg_iso_date()+"'", hdr
   endif
   
return
end
