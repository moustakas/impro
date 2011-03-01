function parse_leda, ledafile, ledapath=ledapath, no_duplicates=no_duplicates
; jm03may19uofa
; we make the assumption that one of the LEDA fields is "NAME"
    
    if n_elements(ledafile) eq 0L then begin
       print, 'Syntax - leda = parse_leda(ledafile,ledapath=)'
       return, -1L
    endif
    
    if n_elements(ledapath) eq 0L then ledapath = cwd()
    
    data = djs_readlines(ledapath+ledafile)
;   data = djs_readlines(ledapath+ledafile,nhead=1,head=head)
    
    dumindx = where(strmatch(data,'!!*',/fold) eq 1B,nhead,comp=dataindx,ncomp=nlines)
    head = data[dataindx[0]]
    data = data[dataindx[1:nlines-1L]]

    nobj = n_elements(data)

    fields = strcompress(strsplit(head,'|',/extract),/remove)
    nfields = n_elements(fields)

;   null = -999.0
    null = ''

    init = create_struct(fields[0],null)
    for i = 1L, nfields-1L do init = create_struct(init,fields[i],null)
    leda = replicate(init,nobj)
    
    for k = 0L, nobj-1L do begin

       line = strsplit(data[k],'|',/extract)
       if strmatch(data[k],'*no cross identification found*',/fold) eq 1B then begin

          name = strn(strmid(data[k],1,strpos(data[k],'<No cross identification')-1))
          leda[k].name = name

;         if n_elements(rem) eq 0L then rem = k else rem = [rem,k]
          
       endif else begin

          for j = 0L, nfields-1L do leda[k].(j) = strcompress(line[j],/remove)
;         for j = 0L, nfields-1L do if strcompress(line[j],/remove) ne '' then leda[k].(j) = line[j]

       endelse

    endfor

;; remove mis-matched objects
;
;    if n_elements(rem) ne 0L then begin
;       
;       good = lindgen(nobj)
;       remove, rem, good
;       leda = leda[good]
;
;    endif

; remove duplicates

    dup = lindgen(n_elements(leda))
    nodup = uniq(leda.name)

    if n_elements(nodup) lt n_elements(dup) then begin

       remove, nodup, dup
       print, 'Duplicate LEDA objects: '
       niceprint, leda[dup].name

       if n_elements(no_duplicates) ne 0L then leda = leda[nodup]

    endif
    
return, leda
end
