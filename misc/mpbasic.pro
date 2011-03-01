;+
; NAME:
;	MPBASIC
; PURPOSE:
;	Provides on-line documentation for IDL topics. The style
;	is a cross between Unix man pages and VMS on-line help. The
;	help is organized in a two level hierarchy --- Level 1 is the
;	global subject, and Level 2 supplies help on subjects within
;	each global subject. If !D.WINDOW is not -1, (window system in use)
;	the mouse is used to prompt for subjects; otherwise, the normal tty
;	interface is used.
;	This routine is used when in Widget IDL and widgets are not available.
; CATEGORY:
;	Help, documentation 
; CALLING SEQUENCE:
;	MPBASIC [, REQUEST]
; INPUTS:
;	REQUEST = A scalar string containing the item for which help is desired.
;	This string can contain one or two (whitespace separated) words.
;	If only one word is included in REQUEST, it is taken as either the 
;       global topic or the program name; if there are two words, the first
;       is taken as the global topic and the second as the program name.   
;       Missing words are prompted for.
; OUTPUTS:
;	Help text is sent to the standard output.
; RESTRICTIONS:
;	The help text is derived from the LaTeX files used to produce
;	the reference manual. However, it is not possible to produce
;	exactly the same output as found in the manual, due to limitations
;	of text-oriented terminals. The text used is therefore considerably
;	abbreviated. Always check the manual if the online help is
;	insufficient. 
;
;	Under VMS, MPBASIC now works (since July 1996) for procedures in text 
;	libraries.    However, it will work faster if the procedures are stored
;	as separate ASCII files in directories.
; MODIFICATION HISTORY:
;	AB, 3, November, 1988						
;	Added ambiguity resolution, ability to handle 
;       multiple levels, and support for mouse.     January, 1989        
;
;       Modified and extended to accept a wider range of help requests,    
;         K.Rhode, STX, June 1990   (Small bug fixed, 9 July 1990)
;       Search the current procedures W. Landsman January, 1991
;	Renamed to MPBASIC for Widget IDL.
;	Modified to match new .HELP file format for IDL 3.1.0 and later. 
;						Joel D. Offenberg, HSTX, 7/1/93
;	Work for VMS text library   W. Landsman    July 1996
;-
function SELECT_TOPIC, SUBJECT, TOPIC_ARRAY, INITIAL
; Given a subject header and an array of topics, returns a string with
; the requested topic (which may or may not be in TOPIC_ARRAY).
; Initial is the index of the initial selection to be highlighted IF  
; a window system menu is used.
on_error,2                      ;Return to caller if an error occurs
xx = fstat(-1)
target = ''
if (!d.name eq 'X' or (!d.window ne -1)) then begin	; Use wmenu
  index = wmenu([SUBJECT, TOPIC_ARRAY,'***CANCEL***'], title=0, initial=initial)
  if (index gt 0) and (index LE N_elements(topic_array)) then $ 
                  target = TOPIC_ARRAY(index-1) else return,''
endif else begin				; Use tty
  if xx.isatty then begin
	  openw, 1, filepath(/TERMINAL), /STREAM, /MORE
	  printf, 1, format = '(/,A,":",/)', SUBJECT
	  printf, 1, TOPIC_ARRAY
	  close, 1
  endif else begin
	  print, format = '(/,A,":",/)', SUBJECT
	  print, TOPIC_ARRAY
  endelse
  print, format='(/,/)'
  read, 'Enter topic for which help is desired: ', target
endelse
target = STRCOMPRESS(STRUPCASE(target),/REMOVE_ALL) ; Upper case & no blanks  
return, target 
end
;
function TOPIC_MATCH, KEY, TOPIC_ARRAY, FOUND, OUTUNIT,EXACT=exact
; Given a string, TOPIC_MATCH returns an array of indices into
; TOPIC_ARRAY that match into FOUND. If there is an exact match
; only its index is returned, otherwise all elements with the same prefix
; match. The number of elements that matched is returned.
; OUTUNIT is the file unit to which output should be directed.
on_error,2                      ;Return to caller if an error occurs
found = [ where(STRTRIM(TOPIC_ARRAY) eq KEY, count) ] ; Match exact string
if (count eq 0) then begin	; No exact match, try to match the prefix
  if keyword_set(EXACT) then return,0
  FOUND = [ where(strpos(TOPIC_ARRAY, KEY) eq 0, count) ]
  if (KEY eq '') then begin
    count = -1
    printf,outunit, !MSG_PREFIX, 'Nothing matching topic "',KEY,'" found.'
    printf,outunit, !MSG_PREFIX, 'Enter "man" for list of available topics.'
  endif
  if ((count gt 1) and (KEY ne '')) then begin
    printf, outunit, format = "(A,'Ambiguous topic ""', A, '"" matches:')",$
     !MSG_PREFIX, KEY
    printf, OUTUNIT, TOPIC_ARRAY(FOUND)
  endif
endif
return, count
end
;
PRO MPBASIC, REQUEST
on_error,2                     
outunit = (inunit = 0)
xx = fstat(-1)
lv1_topic = (lv2_topic = (sv_lv2_topic = (string = ''))) 
m = (no_request = (count1 = (count2 = 0)))
;
if (N_ELEMENTS(REQUEST)) then begin
  temp = size(request)
  if (temp(0) NE 0) then begin
    MSG = 'Argument must be scalar.'
    goto, FATAL
  endif
  if (temp(1) NE 7) then begin
    MSG = 'Argument must be of type string.'
    goto, FATAL
  endif
;
  lv1_topic = STRUPCASE(STRTRIM(STRCOMPRESS(REQUEST), 2))
  sv_lv1_topic = lv1_topic         ; Save the original request for later use
;
;    Parse into one or two strings - level 1 is the global topic, 
;              level 2 the program name
TRY_AGAIN:
  if (((blank_pos = STRPOS(lv1_topic, ' '))) ne -1) then begin
    lv2_topic = STRMID(lv1_topic, blank_pos+1, 10000L)
    lv1_topic = STRMID(lv1_topic, 0, blank_pos)
    if (m eq 0) then sv_lv1_topic = lv1_topic
  endif
endif
;
  ; lv1_files recieves all help files found through !HELP_PATH.
  lv1_dirs = EXPAND_PATH(!HELP_PATH, /ARRAY, COUNT=cnt)
  if (cnt eq 0) then begin
    MSG = 'No online help files found.'
    goto, fatal
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

initial = where(lv1_topics eq 'ROUTINES')
;       If the initial request was simply "man"...
if (m eq 0) and (lv1_topic eq '') then no_request=1
if (lv1_topic eq '') then $
 lv1_topic = SELECT_TOPIC('Help categories', lv1_topics, initial(0)+1)
if lv1_topic eq '' then goto,DONE
if (no_request gt 0) then sv_lv1_topic = lv1_topic
;
if (outunit GT 0) then FREE_LUN, outunit & outunit = 0
if xx.isatty then $
	openw, outunit, filepath(/TERMINAL), /STREAM, /MORE, /GET_LUN $
	else outunit = -1
count = TOPIC_MATCH(lv1_topic,lv1_topics,found,outunit)
;  COUNT tells whether 'lv1_topic' matched any of 'lv1_topics'          
if (m eq 0) then count1 = count      ; if this is the first run-through... 
if (count eq -1) then goto, DONE   ; if no success, then end program
;
;  If the original request was "man", but global topic then requested 
;    was not found, then send a message:

if (no_request gt 0) and (count1 eq 0) then goto,COUNT_ZERO
;
;If the user didn't include a global topic in the initial request, this loop
;adds one, then tries again to find the requested procedure in the help files
M_LOOP:  
if (m le 0) then begin
  if (count eq 0) then begin
    lv2_topic=' '
    lv1_topic='ROUTINES '+STRUPCASE(request)              
    m = m + 1
    if (outunit GT 0) then free_lun, outunit & outunit = 0
    goto,TRY_AGAIN
  endif
endif                                              
if (outunit GT 0) then free_lun, outunit & outunit = 0
lv2_subject = (lv1_files(found))(0)		; Use the first element
;                                        
;If COUNT is STILL zero, a message is sent, and MAN_PROC gives up 
COUNT_ZERO: 
if (count eq 0) then begin    
;
if !Version.os EQ "vms" then begin   
     sep = ','
     dirsep = ''
endif else begin
     sep = ':'
     dirsep = '/'
endelse   
;
temp = !PATH                     ;Get current IDL path of directories
;
;    Loop over each directory in !PATH until procedure name found
;
while temp NE '' do begin   
   dir = gettok(temp,sep)
;
   if strmid(dir,0,1) ne '@' then begin          ;Text Library?
        a = findfile(dir + dirsep +strlowcase(request)+'.pro',COUNT=count)
        if count GE 1 then begin  ;Found by FINDFILE?
       if inunit GT 0 then  close,inunit else get_lun,inunit
        if outunit GT 0 then close,outunit else get_lun,outunit
        OPENR,INUNIT,A(0)
;   Open the output file with the /MORE option
;
if xx.isatty then $
	openw, outunit, filepath(/TERMINAL), /STREAM, /MORE $
	else outunit = -1
        LINE = "STRING"
	WHILE NOT EOF(INUNIT)  DO BEGIN
          READF,INUNIT,LINE
                IF STRMID(LINE,0,2) EQ ';+' THEN BEGIN
                 READF,INUNIT,LINE
                 WHILE NOT EOF(INUNIT) DO BEGIN
      		 READF,INUNIT,LINE
                 IF STRMID(LINE,0,2) EQ ';-' THEN GOTO,DONE
	 PRINTF,OUTUNIT,STRTRIM(STRMID(LINE,1,STRLEN(LINE)))
                ENDWHILE
               ENDIF
	ENDWHILE
        GOTO,DONE
        endif
    endif else begin

         LibName = strmid( Dir, 1, strlen(Dir)-1 )      ;Remove the "@" symbol
         spawn, 'library /list ' + LibName, List
	 lfound = where( list EQ strupcase(request), Nfound)
         if outunit GT 0 then close,outunit else get_lun,outunit
         if Nfound GT 0 then begin
		openw, outunit, filepath(/TERMINAL), /STREAM, /MORE
		spawn,'libr/out=sys$output/extract = ' + request + ' ' + $
			libname, outlist
	 	display = 0
		for i= 0,N_elements(outlist)-1 do begin
			if display then begin
				if strmid(outlist(i),0,2) EQ ';-' then $
					goto,DONE
				printf,outunit,outlist(i)
			endif
	                if strmid(outlist(i),0,2) eq ';+' then display = 1
		endfor
      endif
    endelse
endwhile
endif
if count eq 0 then begin
 if (outunit GT 0) then FREE_LUN, outunit & outunit = 0
  if xx.isatty then $
	openw, outunit, filepath(/TERMINAL), /STREAM, /MORE, /GET_LUN $
	else outunit = -1
  if (sv_lv2_topic ne '') then string = sv_lv1_topic+' '+sv_lv2_topic $
   else string = sv_lv1_topic
  printf,outunit,!MSG_PREFIX, 'Nothing matching topic "',string,'" found.'
  printf,outunit,!MSG_PREFIX, 'Enter "man" for list of available topics.'
  goto,DONE
endif
;
;At this point, a global subject exists - next, process the specific subject
lv2_topics = ''
offset = 0L
if (inunit ne 0) then FREE_LUN, inunit & inunit = 0
;
; Read the specially-formatted .HELP files to look for the requested procedure
  lv2_topics = ''
  offset = 0L
  openr, inunit, lv1_files(found(0)), /GET_LUN
  outunit = 0;
  n = 0L
   tmp=''
  readf,inunit,tmp				; Read first line.
  ; If it's the version tag, parse it.
  version = 1L					; Assume old format
  if (strmid(tmp, 0, 9) eq '%VERSION:') then begin
    reads, tmp, version, format='(9X, I0)'
    readf,inunit,tmp				; Read next line.
  endif
  if (strmid(tmp, 0, 7) eq '%TITLE:') then readf, inunit, tmp   ; Skip title
  n = long(tmp)					;# of records
  if (version ne 1) then begin
    ; Version 2 format has the number of characters used by all the
    ; subtopics on the next line. We don't use it, but have to read it
    readf,inunit,tmp				; Read next line.
  endif

;
;	Search the beginning of the .HELP file for a number.  Disregard
;	all lines before the first number (added for compatibility with
;	IDL 3.1.0's new .HELP file format).
;	


;WHILE (n eq 0) do BEGIN
;	dummy = ""
;	readf, inunit, dummy 
;	dd = byte(strcompress(dummy,/REMOVE_ALL))
;
;	;Is the first character a number (between '0' and '9')?
;	IF (dd(0) GE 48 and dd(0) LE 57 ) THEN $
;		n = fix(dummy)
;endWHILE

;The following line is no longer needed.
;readf,inunit,n			        ;Read # of records

lv2_topics = strarr(n)			;Make names
readf,inunit,lv2_topics			;Read entire string to inunit
if (version EQ 1) then begin
	offsets = long(strmid(lv2_topics, 15, 30))	;Extract starting bytes
	lv2_topics = strmid(lv2_topics,0,15)		;Isolate names
  endif else begin
    offsets = lonarr(n)
    for i = 0, n-1 do begin
      tmp = lv2_topics(i)
      colon = strpos(tmp, ':') + 1		; Find delimiter
      offsets(i) = long(strmid(tmp, 0, colon))
      lv2_topics(i) = strmid(tmp, colon, 10000000)
    endfor
  endelse
tmp = fstat(inunit)        ; Determine the base of the help text in .HELP files
text_base = tmp.cur_ptr
if text_base eq tmp.size then category = 1  else category =0
;
; If no level 2 topic has been supplied, prompt for one
if lv2_topic eq '' then $
 lv2_topic = SELECT_TOPIC(STRUPCASE(lv2_subject), lv2_topics, 1)
if lv2_topic eq '' then goto,DONE
;
if (m eq 0) then sv_lv2_topic = lv2_topic 
if (outunit GT 0) then FREE_LUN, outunit & outunit = 0
if xx.isatty then $
	openw, outunit, filepath(/TERMINAL), /STREAM, /MORE, /GET_LUN $
	else outunit = -1
;
; If count is still zero, and all the ROUTINES have been searched, then quit
if (((count=TOPIC_MATCH(lv2_topic,lv2_topics,found,outunit,/EXACT))) eq 0) $
 and (m eq 1) then goto,COUNT_ZERO
count2 = count
;
;If the global topic was found, but the procedure requested doesn't exist, quit
if (count2 eq 0) and (count1 eq 1) then goto,COUNT_ZERO
;
; If all of the ROUTINES been searched, go back to M_LOOP
if (count eq 0) and (m eq 0) then goto,M_LOOP 
;
; Print the documentation if the search was successful
str = ''
if (outunit GT 0) then FREE_LUN, outunit & outunit = 0
if category then begin
   if (inunit NE 0) then FREE_LUN,inunit & inunit = 0      ;Corrected Sep 91
   count = 0 & m= -1
   request = lv2_topic
   goto,M_LOop
endif
if xx.isatty then $
	openw, outunit, filepath(/TERMINAL), /STREAM, /MORE, /GET_LUN $
	else outunit = -1
for i=0,(count-1) do begin
  index = found(i)
  if (count gt 1) then printf, outunit, lv2_topics(index), $
   format='("***************",/,A,/,"***************")'
  POINT_LUN, inunit, text_base + offsets(index)
  readf, inunit, str		; Skip the ";+"
  readf, inunit, str
  while (STRTRIM(str) NE ";-") do begin
    printf, outunit, str, ' '
    readf, inunit, str
  endwhile
endfor
goto,DONE
;
FATAL:		; The string MSG must be already set
message, MSG   
;
DONE:
if (outunit GT 0) then FREE_LUN, outunit & outunit = 0
if (inunit ne 0) then FREE_LUN, inunit & inunit = 0
end
