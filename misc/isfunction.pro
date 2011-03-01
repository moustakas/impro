function isfunction,proname, outnames, numline
;+
; NAME:
;	ISFUNCTION()
; PURPOSE:
;	Determine whether the IDL program(s) in a file are procedures or 
;	functions.    Needed because the intrinisc RESOLVE_ROUTINE and 
;	ROUTINE_INFO() procedures require the user to know beforehand whether 
;	to supply the /IS_FUNCTION or /FUNCTION keywords.
;
; CALLING SEQUENCE:
;	status = ISFUNCTION( filename, [ outnames, numlines]
; INPUT:
;	filename = scalar string giving complete specification if file name
;		(include .pro extension)
;
; OUTPUT:
;	status - integer vector with number of elements equal to the number 
;	of routines in the file.    Each status value consists of 0 or 1
;	 1 - routine is an IDL function
;	 0 - routine is an IDL procedure
;	 If no valid IDL functions or procedures are found in the file, then
;		ISFUNCTION() returns a scalar value of -1 
;
; OPTIONAL OUTPUTS:
;	outnames - vector string, giving name of each IDL procedure or function
;		in the file
;	numlines - integer vector, giving the number of lines in each IDL
;		procedure or function in the file
; PROCEDURE CALLS:
;	FDECOMP
; REVISION HISTORY:
;	Written, W. Landsman                  June, 1995

 openr,lun,proname,/get_lun
 line = ''
 FDECOMP, proname, disk, dir, name, ext
 pname = strtrim(strlowcase(name),2)
 outnames = ''
 status = -1
 numline = 0
 nprog = 0
 step = 0
 while not eof(lun) do begin
	readf,lun,line
	comment = strpos(line,';')
	if comment EQ 0 then goto, Next_line
	if comment GT 0 then line = strmid(line,0,comment-1)
	tabpos = strpos(line,string(9b))    ;Remove all tabs
	while tabpos ne -1 do begin
		strput,line,' ',tabpos
		tabpos = strpos(line,string(9b) )
	endwhile
	line = ' ' + strlowcase(line)
		
	pos = strpos(line, ' pro ')
	if ( pos GE 0 ) then begin
		comma = strpos(line,',')
		if comma GT 0 then line = strmid(line,0,comma)
		procname = strmid(line,pos+5,80)
		outnames = [outnames,procname]
		numline = [numline,step]
		status = [status,0]
                step = 0
       endif else begin
		
        pos = strpos(line,' function ')
	if ( pos GE 0 ) then begin
		comma = strpos(line,',')
		if comma GT 0 then line = strmid(line,0,comma)
		procname = strmid(line,pos+10,80)
		outnames = [outnames,procname]
		numline = [numline, step]
		status = [status,1]
		step = 0
	endif
	endelse 
NEXT_LINE:
	step  = step +1

 endwhile
 free_lun,lun

; Any comments at the top of the file are included in the line count for the
; first procedure.

 if N_elements(status) GT 1 then begin
      numline = [numline,step]
      numline[2] = numline[1] + numline[2]    
      nprog = N_elements(status) - 1
      If N_elements(status) GT 1 then begin 
	outnames = strtrim( strlowcase(outnames[1:*]),2)
 	status = status[1:*]
	numline = numline[2:*]
      endif
 endif
 
 return,status
 end
