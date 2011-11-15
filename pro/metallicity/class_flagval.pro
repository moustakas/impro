;+
; NAME:
;   class_flagval
;
; PURPOSE:
;   Return bitmask values corresponding to labels in support of
;   ICLASSIFICATION. 
;
; CALLING SEQUENCE:
;   value = class_flagval(flagprefix, label)
;
; INPUTS:
;   flagprefix - Flag name (scalar string)
;   label      - String name(s) corresponding to each non-zero bit in FLAGVALUE.
;
; OPTIONAL KEYWORDS:
;   /check     - checking if a good flag name, so don't abort if not
;
; OUTPUTS:
;   value      - Signed long with any number of its bits set.
;
; COMMENTS:
;   This function is the inverse of CLASS_FLAGNAME().
;
; PROCEDURES CALLED:
;   splog
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IMPRO_DIR/etc/class_maskbits.par
;
; REVISION HISTORY:
;   02-Apr-2002 Written by D. Schlegel, Princeton.
;   29-Jul-2008 imported into IMPRO, Moustakas
;-

function class_flagval, flagprefix, inlabel, check=check

; Declare a common block so that the mask names are remembered between calls.

;   common com_im_maskbits, maskbits

    if (n_params() NE 2 OR n_elements(flagprefix) NE 1) then begin
       print, 'Syntax - value = class_flagval(flagprefix, label)'
       return, ''
    endif

; Read the parameter file the 1st time this function is called.
; (After that, store this info in a common block.)

    maskfile = filepath('class_maskbits.par',root_dir=getenv('IMPRO_DIR'),$
      subdirectory='etc')
    if (file_test(maskfile,/regular) eq 0L) then $
      message, 'File with mask bits not found'
    yanny_read, maskfile, pdat
    maskbits = *pdat[0]
    yanny_free, pdat

;   if (NOT keyword_set(maskbits)) then begin
;      maskfile = sg1120_path(/sex)+'im_maskbits.par'
;      if (NOT keyword_set(maskfile)) then $
;        message, 'File with mask bits not found'
;      yanny_read, maskfile, pdat
;      maskbits = *pdat[0]
;      yanny_free, pdat
;   endif

; Generate a list of all non-blank labels as a string array

    flagvalue = 0L

    alllabel = strsplit(inlabel[0], /extract)
    for i=1, n_elements(inlabel)-1 do $
      alllabel = [alllabel, strsplit(inlabel[i], /extract)]
    ilabel = where(alllabel NE '', nlabel)
    if (nlabel EQ 0) then return, flagvalue
    alllabel = alllabel[ilabel]

; Find the match for each label, and add its value to the output

    for ilabel=0, nlabel-1 do begin
       imatch = where(strupcase(flagprefix[0]) EQ maskbits.flag $
         AND strupcase(alllabel[ilabel]) EQ strupcase(maskbits.label), ct)
       if (ct NE 1) then begin
          if(NOT keyword_set(check)) then begin
             message, 'ABORT: Unknown bit label ' + $
               strupcase(alllabel[ilabel]) $
               + ' for flag ' + strupcase(flagprefix)
          endif else begin
             return, 0
          endelse
       endif

       flagvalue = flagvalue + 2L^(maskbits[imatch[0]].bit)
    endfor

return, flagvalue
end
