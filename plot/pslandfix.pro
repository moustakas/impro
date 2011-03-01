pro pslandfix,filename, xw=xw, yw=yw
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

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: pslandfix,filename, xw=xw, yw=yw'
      return
  ENDIF 
  
  IF n_elements(yw) EQ 0 THEN yw=560
  IF n_elements(xw) EQ 0 THEN xw=40
  stop
  putstr = ntostr(long(yw))+' '+ntostr(long(xw))
  
  openr,1,filename
  lin=''
  replinz=strarr(100)
  linno=lonarr(100)
  ctr=0 & linctr=0L
  
  IF !version.release LT 5.4 THEN BEGIN
      search_string = 'save $IDL_DICT begin'
      sr = [0,20]
  ENDIF ELSE BEGIN 
      search_string = '$IDL_DICT begin'
      sr = [0,15]
  ENDELSE 

  print,'Fixing Landscape for file: '+filename
  while not EOF(1) do begin
      readf,1,lin
      if (strmid(lin,sr[0],sr[1]) eq search_string) then begin
          print,lin
          
          strput,lin,putstr,sr[1]+1  ;mine 2
;      strput,lin,'576 20',21    ;Deutch
;          strput,lin,' 90',58
          strput,lin,' 90',sr[1]+1+37
          print,lin
          print,'---'
          replinz(ctr)=lin
          linno(ctr)=linctr
          ctr=ctr+1
      endif
      linctr=linctr+1
  endwhile
  close,1
  
  
  openu,1,filename
  lin=''
  ctr=0 & linctr=0L
  
  while not EOF(1) do begin
      if (linctr ne linno(ctr)) then begin
          readf,1,lin
      endif else begin
          printf,1,replinz(ctr)
          ctr=ctr+1
      endelse
      linctr=linctr+1
  endwhile
  close,1
  
  return
end
