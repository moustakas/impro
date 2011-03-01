;+
; NAME: ARM_NAMEPARSE
;       
; PURPOSE: parse a file name for root, extension and path
;
; CATEGORY: file management
;
; CALLING SEQUENCE: ARM_NAMEPARSE, files, [root=, extension=, path=]
;
; INPUTS: files - array of one or more file names (can include path)
; 
; OPTIONAL OUTPUTS: 
;    root      - root name of file(s)
;    extension - file extension(s)
;    path      - file path(s)
;
; EXAMPLE: 
;    IDL> fname = '/home/psfiles/idl.ps'
;    IDL> arm_nameparse, fname, root=root, extension=ext, path=path
;
;    root -> idl, ext -> ps, path -> /home/psfiles/
;
; ROUTINES: ARM_ERROR
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., September 2003
;    now accepts array of files, ARM, October 16, 2003
;-

pro arm_nameparse, files, root=root, extension=extension, path=path

; initialize the output variables as null strings

    n = N_ELEMENTS(files)
    root      = STRARR(n)
    extension = STRARR(n)
    path      = STRARR(n)

    for i=0,n-1 do begin

       file = files[i]

       parts = STRSPLIT(file, '/', /regex, /extract)
       nparts = N_ELEMENTS(parts)

; if no backslash, no path

       if nparts eq 1 then fname = file $
       else begin
          
; if backslash, identify path
          
          path[i] = STRJOIN([parts[0:nparts-2],''], '/')
          if STRPOS(file, '/') eq 0 then path[i] = '/'+path[i]
          fname = parts[nparts-1]
       endelse
       
       lastslash = STRPOS(file, '/', /reverse_search)
       if lastslash ne STRLEN(file)-1 then begin
          
; if no trailing backslash, identify root/extension
          
          parts = STRSPLIT(fname, '.', /extract)
          nparts = N_ELEMENTS(parts)
          
          if nparts gt 1 then begin
             
; if root and extension separated by period
             
             extension[i] = parts[nparts-1]
             root[i] = STRJOIN(parts[0:nparts-2], '.')
             
          endif else $
            
; if no period, no extension
          
          if not STRMATCH(fname, '*.*') then root[i] = parts $
            
; if period but no root
          
          else if STRMATCH(fname, '.*') then extension[i] = parts $
            
; if period but no extension
          
          else root[i] = parts
          
       endif $
         
; if no root or extension, reappend trailing backslash
       
       else if lastslash ne -1 then path[i] = path[i] + parts[nparts-1] + '/'
       
    endfor
    
    return
    
 end
