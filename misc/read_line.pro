
pro read_line, fp, line, comment_string=cstr, notrim=notrim, help=help

;+
; NAME:
;    read_line 
;
; PURPOSE:
;    Simple procedure to read a line from a file, ignoring comments and
;    blank lines.  Does nothing if the end of file has already beem
;    reached
;
; CATEGORY:
;    Input/Output Routines
;
; CALLING SEQUENCE:
;    read_line, fp, line[, comment_str, notrim, help]
; 
; INPUTS:
;    fp             Unit number of the input file
;
; OPTIONAL INPUTS:
;	
; KEYWORD PARAMETERS:
;    comment_string A string containing all allowed symbols which may
;                   denote a comment if found at the beginning of a
;                   line.  Defaults to '#'
;
;    notrim         Set this keyword to disable string compression
;
; OUTPUTS:
;    line           The string read from the file   
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;    Strips line of leading/trailing whitespaces, and compresses other
;    whitespaces into single ones
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;    read_line, 'foo.dat', line
;
; MODIFICATION HISTORY:
;    1996aug06  JEB  Written
;    1996sep09  JEB  Modified to allow no string trimming
;-

if n_params() lt 2 or keyword_set(help) then begin
    doc_library, 'read_line'
    retall
endif

if not eof(fp) then begin
    if not keyword_set(cstr) then cstr = '#'
    line = ' '

    readf, fp, line
    line1 = strcompress(strtrim(line, 2))

; --Skip blanks and comment lines
    while (strlen(line1) eq 0 or strpos(cstr, strmid(line1,0,1)) ne -1) $
      and not eof(fp) do begin
        readf, fp, line
        line1 = strcompress(strtrim(line, 2))
    endwhile
    if not keyword_set( notrim ) then line = line1
endif 

end
