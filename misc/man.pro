;+
; NAME:
;	MAN
; PURPOSE:
;	Provides online documentation for IDL topics. If the current
;	graphics device supports widgets, a graphical user interface
;	is used. Otherwise, a more basic version which is a cross
;	between Unix man pages and VMS online help is used. 
;	This routine can emulate the V3.x online help system in V4.x and later.
; CATEGORY:
;	Help, documentation .
; CALLING SEQUENCE:
;	MAN [, REQUEST]
; INPUTS:
;	REQUEST = A scalar string containing the procedure name or topic
;		  for which help is desired.    MAN first tries to interpret
;		  the string as a topic, and if not found then searches the
;		  !PATH for the procedure.
; OUTPUTS:
;	The widget version uses a text widget to display the help
;	text. The basic version sends help text to the standard output.
;
; EXAMPLES:                                                             
;                             
;	IDL> man,'rebin'   ;Display help on intrinisic IDL REBIN function
;	IDL> man,'hrebin'  ;Display help on Astronomy Library HREBIN procedure
;	IDL> man,'astron'  ;Display all astronomy library procedures
; NOTES:
;	To install the MAN procedure
;
;	1.  Place the procedures in the /contrib/landsman directory in one's 
;	!PATH.   
;
;	2.  Define the !HELP_PATH sytem prior to entering IDL, to point to a 
;	directoryto contain the help files.    This can be done either by 
;	setting the environment variable IDL_HELP_PATH or by assigning 
;	!HELP_PATH in a startup file,
;
;	3.  Place any IDL help files (such ASTRON.HELP  from the Astronomy 
;	Library) in the !HELP_PATH directory.    Also place the file 
;	ROUTINES.HELP from IDL V3.6 into this directory.   This information 
;	in the ROUTINES.HELP file will be somewhat obsolete, but it is the 
;	only way to get online help for intrinsic procedures.
;
; RESTRICTIONS:
;	The help text for intrinsic is derived from the LaTeX files used to 
;	produce the V3.6 reference manual. However, it is not possible to 
;	produce exactly the same output as found in the manual due to 
;	limitations of text oriented terminals. The text used is therefore 
;	considerably abbreviated. Always check the manual if the online help is
;	insufficient. 
;
;	MAN cannot locate procedures in a VMS text library
; MODIFICATION HISTORY:
;	Adapted from MAN_PROC, W. Landsman             February, 1996
;-
;

PRO MAN, REQUEST

  ; lv1_files recieves all help files found through !HELP_PATH.
  lv1_dirs = EXPAND_PATH(!HELP_PATH, /ARRAY, COUNT=cnt)
  if (cnt eq 0) then begin
     message,'No online help files found.',/CONT
     return
  endif
  for i = 0, cnt-1 do begin
    tmp = STRLOWCASE(findfile(filepath('*.help', root_dir=lv1_dirs(i))))
    if (i eq 0) then lv1_files = TEMPORARY(tmp) $
    else lv1_files=[lv1_files, TEMPORARY(tmp)]
  endfor

  ; lv1_topics gets uppercase version of just the names.
  lv1_topics = STRUPCASE(lv1_files)
  stlen = strlen(lv1_topics)
  if !version.os eq 'windows' then begin 
    tail = STRPOS(lv1_topics, '.HEL')
  endif else if !version.os eq 'vms' then begin 
    tail = STRPOS(lv1_topics, '.HELP;')
  endif else $
    tail = STRPOS(lv1_topics, '.HELP')
  n = n_elements(lv1_topics)
  for i = 0, n-1 do $
	lv1_topics(i) = strmid(lv1_topics(i), 0, tail(i))
  for i = 0, n-1 do begin	; Strip path part off lv1_topics
    case !version.os of
      'vms': begin
           j = STRPOS(lv1_topics(i), ']')
           while (j ne -1) do begin
             lv1_topics(i) = strmid(lv1_topics(i), j+1, 32767)
             j = STRPOS(lv1_topics(i), ']')
           endwhile
      end
      'windows': begin
        j = STRPOS(lv1_topics(i), '\')
        while (j ne -1) do begin
  	  lv1_topics(i) = strmid(lv1_topics(i), j+1, 32767)
          j = STRPOS(lv1_topics(i), '\')
        endwhile
      end
      'MacOS': begin
        j = STRPOS(lv1_topics(i), ':')
        while (j ne -1) do begin
  	  lv1_topics(i) = strmid(lv1_topics(i), j+1, 32767)
          j = STRPOS(lv1_topics(i), ':')
        endwhile
      end
      else:  begin      ; Unix otherwise
        j = STRPOS(lv1_topics(i), '/')
        while (j ne -1) do begin
  	  lv1_topics(i) = strmid(lv1_topics(i), j+1, 32767)
          j = STRPOS(lv1_topics(i), '/')
        endwhile
      end
    endcase
  endfor

  ; Sort the topics into alphabetical order.
  tmp = sort(lv1_topics)
  lv1_files = lv1_files(tmp)
  lv1_topics = lv1_topics(tmp)

;Determine if a request is present.  If a request is present, determine if it 
;is a Level 1 topic (i.e. .HELP file) or a Level 2 topic (i.e. IDL routine).
;If the latter, call mpBasic, not widgets_olh.

    if N_elements(Request) EQ 0 then Request = ''
    SZREQ = size(Request)
    lnreq = strlen(Request)

    if SZREQ(0) ne 0 then begin
	message,/inf,'Request must be a scalar'
	return
    endIF ELSE IF (szreq(1) ne 7) AND (szreq(1) ne 0) then BEGIN
	message,/inf,'Request must be a string'
	return
    endIF ELSE IF (lnreq eq 0) then BEGIN		;No request present
	  IF ((!D.FLAGS and 65536) eq 65536) then $	;Widgets present?
		WIDGET_OLH, REQUEST $
	  ELSE MPBASIC, REQUEST
	return
    endIF ELSE BEGIN			;Request present. 
 	;Is there a space?  If yes, then the two-word request goes to MPBASIC.
	IF ((strpos(request,string(32b)) gt 0) and $
	    (strpos(request,string(32b)) lt strlen(request)-1)) THEN BEGIN
		MPBASIC, REQUEST
		return
	endIF
	;Is the word a Level 1 topic?  
	IF (total(where(lv1_topics eq strupcase(strcompress(request,$
	/REMOVE_ALL)))) GE 0) then BEGIN
		IF ((!D.FLAGS and 65536) eq 65536) then $     ;Widgets present?
			WIDGET_OLH, REQUEST $
		ELSE MPBASIC, REQUEST
		return
	endIF ELSE $
		MPBASIC, REQUEST
    endELSE
end
