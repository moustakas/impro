;+
; NAME: arm_maketable
;       
; CATEGORY: miscellaneous
;
; PURPOSE: generate table from unformatted, delineated lines of text 
;
; CALLING SEQUENCE: ARM_MAKETABLE, infile, [outfile=, delimiter=, 
;                     new_delimiter=, comment=, ignore=, skipline=,
;                     /skipblanks]
;
; INPUTS: 
;   infile - input file name (string)
;       
; OPTIONAL INPUTS:
;   outfile       - name of output file (string), 
;                   default=arm_maketable.out
;   delimiter     - character used to separate column entries in INPUT 
;   new_delimiter - string used to separate column entries in OUTPUT
;   comment       - string character used to comment lines in INPUT
;   ignore        - disregard INPUT lines beginning with this character
;   skipline      - number of initial lines to skip in INPUT
;
; KEYWORDS:
;   skipblanks - ignore blank lines
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; COMMENTS: 
;
; BUG REPORT: Please report any bugs to Andrew R. Marble.
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs.
;-

pro ARM_MAKETABLE, infile, outfile=outfile, delimiter=delimiter, $
       alignment=alignment, comment=comment, ignore=ignore, $
       skipblanks=skipblanks, skipline=skipline, new_delimiter=new_delimiter

; defaults and error checking

    if FILE_SEARCH(infile) eq '' then begin
       MESSAGE, 'THE FILE '+infile+' COULD NOT BE LOCATED.', /continue
       return
    endif

    if N_ELEMENTS(outfile) eq 0 then outfile = 'arm_maketable.out'
    if FILE_SEARCH(outfile) ne '' then begin
       MESSAGE, 'THE FILE '+outfile+' ALREADY EXISTS. WRITE-OVER? (y/n).', /continue
       answer = STRUPCASE(GET_KBRD(1))
       while TOTAL(STRMATCH(['Y','N'],answer)) eq 0 do answer = STRUPCASE(GET_KBRD(1))
       if answer eq 'N' then begin
          MESSAGE, 'ARM_MAKETABLE ABORTED.'
          return
       endif
    endif
    
    if N_ELEMENTS(skipline) eq 0 then skipline = 0
    if N_ELEMENTS(delimiter) eq 0 then delimiter = ' '
    if N_ELEMENTS(comment) eq 0 then comment = 'dummytext=nocomment'
    if N_ELEMENTS(new_delimiter) eq 0 then new_delimiter = ''
    if new_delimiter ne '' then new_delimiter = ' ' + new_delimiter + ' '

; read in table rows

    line = ''
    first = 1L

    OPENR, lun, infile, /get_lun
    counter = 0L
    while not EOF(lun) do begin
       READF, lun, line
       if counter ge skipline then begin
          if first then begin
             lines = line
             first = 0L
          endif else lines = [lines, line]
       endif
       counter = counter + 1L
    endwhile
    CLOSE, lun
    FREE_LUN, lun

; discard ignored lines

    bad = BYTARR(N_ELEMENTS(lines))
    for i = 0, N_ELEMENTS(ignore)-1 do begin
       bad = bad + STRMATCH(STRMID(STRTRIM(lines,2),0,1),ignore[i]) < 1
       if KEYWORD_SET(skipblanks) then bad = $
         bad + STRMATCH(lines,'') < 1 
    endfor
    if TOTAL(bad) gt 0 then REMOVE, WHERE(bad), lines

    commented = BYTARR(N_ELEMENTS(lines))

    first = 1L

; separate column values in uncommented lines

    for i = 0, N_ELEMENTS(lines)-1 do begin

       if TOTAL(STRMATCH(comment,STRMID(STRTRIM(lines[i],2),0,1))) eq 0 and $
         lines[i] ne '' then begin
            
          fields = STRTRIM(STRSPLIT(lines[i],delimiter,/extract,/regex),2)
          if first then begin
             cols = fields 
             first = 0L
          endif else cols = [[cols],[fields]]

       endif else commented[i] = 1
             
    endfor

; make columns same length

    ncols = N_ELEMENTS(cols[*,0])
    for i = 0, ncols-1 do cols[i,*] = ARM_STRFMT(cols[i,*])

; output table

    OPENW, lun, outfile, /get_lun

    cols[0:ncols-2,*] = cols[0:ncols-2,*] + new_delimiter
    
    counter = 0L
    fmt = '('+STRJOIN(REPLICATE('A',ncols),',')+')'
    for i = 0, N_ELEMENTS(lines)-1 do begin
       if commented[i] then PRINTF, lun, f='(A)', lines[i] else begin
          PRINTF, lun, f=fmt, cols[*,counter]
          counter = counter + 1
       endelse
    endfor
    CLOSE, lun
    FREE_LUN, lun

end
