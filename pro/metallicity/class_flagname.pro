;+
; NAME:
;   class_flagname
;
; PURPOSE:
;   Return bitmask labels corresponding to bit numbers in support of
;   ICLASSIFICATION. 
;
; CALLING SEQUENCE:
;   label = class_flagname(flagprefix, flagvalue, [ /concat, /silent ] )
;
; INPUTS:
;   flagprefix - Flag name (scalar string).  The following are supported:
;                PIXMASK
;   flagvalue  - Signed long with any number of its bits set.
;
; OPTIONAL KEYWORDS:
;   concat     - If set, then concatenate all of the output labels in
;                LABEL into a single whitespace-separated string.
;   silent     - If set, then don't print a warning when there is no bit label
;                corresponding to one of the bit values.
;
; OUTPUTS:
;   label      - String name(s) corresponding to each non-zero bit in FLAGVALUE.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This function is the inverse of CLASS_FLAGVAL().
;
; EXAMPLES:
;
; BUGS:
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
;   01-Apr-2002 Written by D. Schlegel, Princeton.
;   29-Jul-2008 imported into IMPRO, Moustakas
;-
;------------------------------------------------------------------------------
function class_flagname, flagprefix, flagvalue, concat=concat, silent=silent

   ; Declare a common block so that the mask names are remembered between calls.
;  common com_im_maskbits, maskbits

   if (n_params() NE 2 OR n_elements(flagprefix) NE 1) then begin
      print, 'Syntax - label = class_flagname(flagprefix, flagvalue, [ /concat ] )'
      return, ''
   endif

   ;----------
   ; Read the parameter file the 1st time this function is called.
   ; (After that, store this info in a common block.)

    maskfile = filepath('class_maskbits.par',root_dir=getenv('IMPRO_DIR'),$
      subdirectory='etc')
    if (file_test(maskfile,/regular) eq 0L) then $
      message, 'File with mask bits not found'
    yanny_read, maskfile, pdat
    maskbits = *pdat[0]
    yanny_free, pdat

   ;----------
   ; Find the match for each non-zero bit.

   indx = where(djs_int2bin(flagvalue), nret)
   if (indx[0] EQ -1) then begin
      retval = ''
   endif else begin
      retval = strarr(nret)
      for iret=0, nret-1 do begin
         j = where(strupcase(flagprefix[0]) EQ maskbits.flag $
          AND indx[iret] EQ maskbits.bit)
         if (j[0] NE -1) then retval[iret] = maskbits[j].label $
          else if (NOT keyword_set(silent)) then $
           splog, 'MESSAGE: Unknown bit ', indx[iret], $
           ' for flag ' + strupcase(flagprefix)
      endfor
   endelse

   ;----------
   ; If /CONCAT is set, then concatenate all of the output strings
   ; into a single string separted only by whitespace.

   if (keyword_set(concat)) then begin
      for i=1, nret-1 do $
       retval[0] = retval[0] + ' ' + retval[i]
      retval = retval[0]
   endif

   return, retval
end
;------------------------------------------------------------------------------
