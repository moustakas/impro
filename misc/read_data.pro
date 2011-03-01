
pro read_data, file, x, y, $
               columns=cols, ncolumns=ncols, $
               row_start=row1, nrows=nrows, $
               comment_string=cstr, silent=silent, help=help

;+
; NAME:
;    read_data
;
; PURPOSE:
;    Reads columnar data from an ASCII file
;
; CATEGORY:
;    Input/Output Routines
;
; CALLING SEQUENCE:
;    read_data, file, x(, y[, columns, ncolumns, row_start, nrows, 
;               comment_str, silent])
; 
; INPUTS:
;    file           Input file name
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;    columns        Array containing column numbers to read.  The first
;                   column in the file is column one
;
;    ncolumns       If the columns keyword is not set, this specifies the
;                   number of columns to read, starting from the left.
;                   Otherwise it has no effect
;
;    row_start      Start reading from this row of data (not including
;                   comment/blank lines).  Defaults to the first row
;
;    nrows          Read this many rows starting from row_start.  Default
;                   is all rows in the file
;
;    comment_string A string containing all allowed symbols which may
;                   denote a comment if found at the beginning of a
;                   line.  Defaults to '#'
;
;    silent         If set, does not output any diagnostic information
;
; OUTPUTS:
;    x               Output array of data.  If y is given as an argument, x
;                    will be a one-dimensional array containing the first
;                    column read.  Otherwise, the ith column is stored
;                    in x(*, i-1)
;
; OPTIONAL OUTPUTS:
;    y               If y is given as an argument, the first column read is 
;                    stored in x(*) and subsequent columns i>1 are stored 
;                    in y(*, i-2)
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;    Overwites the arguments x and y
;
; RESTRICTIONS:
;    Assumes that columns are delimited by whitespaces in the input file.
;    The procedure will crash if one of the requested columns for reading
;    contains non-numeric gobbledygook (outside of comment lines)
;
; PROCEDURE:
;    Requires the procedures get_words and read_line.  This is a 
;    generalised version of read2.pro by M. Liu, with many added 
;    keywords.  In addition, this routine runs somewhat faster in the
;    default case where all columns are read in order, because the
;    inner loops have been cleaned up.  Columns should be specified
;    only when many of the columns in the file are not to be read, or
;    if any column contains non-numeric text, because the procedure is
;    somewhat slower when columns are specified
;
; EXAMPLES:
;    read_data, 'foo.dat', x, y, col=[2,4,5], comm='!'
;    read_data, 'foo.dat', x, ncol=3, /silent
;    read_data, 'foo.dat', x, y
;    read_data, 'foo.dat', x, col=[2,1], row_start=100, nrows=50
;
; MODIFICATION HISTORY:
;    1996aug06  JEB  Written
;    1996sep13  JEB  Bug fix
;    1999nov22  JEB  Use long int for row index
;-

On_error, 2

if n_params() lt 2 or keyword_set(help) then begin
    doc_library, 'read_data'
    retall
endif
nvar = n_params() - 1

; --Default keyword parameters
if not keyword_set(row1) then row1 = 1L
if not keyword_set(cstr) then cstr = '#'

; --Input file
openr, ifp, file, /get_lun

; --Skip to row specified by row_start
for irow = 1, row1-1 do read_line, ifp, line, comm=cstr

; --Read first line
nr = row1
line = ' '
read_line, ifp, line, comm=cstr

; --Count number of data columns in first line
get_words, line, words
nc = n_elements(words)

; --Check that data are not too small
warnstr = 'data have too few columns for array in readdata.'

; --Must be able to read all the columns requested
if keyword_set(cols) then begin
    ncols = n_elements(cols) 
    if nc lt max(cols) then begin
        print, 'Error: ' + warnstr
        retall
    endif
endif else if keyword_set(ncols) then begin
    if nc lt ncols then begin
        print, 'Error: ' + warnstr
        retall
    endif 

; --Default: read all columns found in the first line
endif else ncols = nc

; --If y is defined, need at least 2 columns
if nvar gt 1 and ncols lt 2 then begin
    print, 'Error: ' + warnstr
    retall
endif

; --Count remaining non-comment, non-blank rows
while not eof(ifp) do begin
    read_line, ifp, line, comm=cstr
    if strlen(line) ne 0 then nr = nr + 1
endwhile

; --Not enough rows to read in all requested
if keyword_set(nrows) then begin
    if (row1+nrows-1 gt nr and not keyword_set(silent)) then begin
        nrows = nr - row1 + 1
        print, format = '(" Warning: reading only ", i0, " rows!")', nrows
    endif
endif else nrows = nr - row1 + 1

free_lun, ifp

if not keyword_set(silent) then $
  print, format = '("...Reading ", i0, " row(s) and ", i0, " column(s)...")', $
  nrows, ncols

; --Initialise arrays
if nvar eq 1 then begin
    x = fltarr(nrows, ncols)
    xc = fltarr(ncols)
endif else begin
    x = fltarr(nrows)
    y = fltarr(nrows, ncols-1)
    yc = fltarr(ncols-1)
endelse

; --Now read in the data
openr, ifp, file, /get_lun

; --Skip to row specified by row_start
for irow = 1, row1-1 do read_line, ifp, line, comm=cstr

; --Use columns defined by keyword parameter
if keyword_set(cols) then begin

; --Store in x array
    if nvar eq 1 then begin
        irow = 0L
        read_line, ifp, line, comm=cstr
        while irow lt nrows do begin
            get_words, line, words, nwords=max(cols)
            for icol = 0, ncols-1 do begin
                iw = cols(icol) - 1
                reads, words(iw), xx
                x(irow, icol) = xx
            endfor
            irow = irow + 1
            read_line, ifp, line, comm=cstr
        endwhile

; --Store in x and y arrays
    endif else begin
        irow = 0L
        read_line, ifp, line, comm=cstr
        while irow lt nrows do begin
            get_words, line, words, nwords=max(cols)
            iw = cols(0) - 1
            reads, words(iw), xx
            x(irow) = xx
            for icol = 1, ncols-1 do begin
                iw = cols(icol) - 1
                reads, words(iw), xx
                y(irow, icol-1) = xx
            endfor
            irow = irow + 1
            read_line, ifp, line, comm=cstr
        endwhile
    endelse

; --Take columns in order
endif else begin

; --Store in x array
    if nvar eq 1 then begin
        irow = 0L
        read_line, ifp, line, comm=cstr
        while irow lt nrows do begin
            reads, line, xc
            x(irow, *) = xc
            irow = irow + 1
            read_line, ifp, line, comm=cstr
        endwhile

; --Store in x and y arrays
    endif else begin
        irow = 0L
        read_line, ifp, line, comm=cstr
        while irow lt nrows do begin
            reads, line, xc, yc
            x(irow) = xc
            y(irow, *) = yc
            irow = irow + 1
            read_line, ifp, line, comm=cstr
        endwhile
    endelse

endelse

free_lun, ifp
 
end
