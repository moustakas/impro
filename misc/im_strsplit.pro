function im_strsplit, str, key, _extra=extra
; jm09may21nyu - simple wrapper on strsplit() that works on vectors;
; for example, if STR=['2006-Jan-10','2009-Feb-10'] and KEY='-' then: 
;   IDL> split = im_strsplit(str,key,/extract)
;   IDL> help, split  
;   SPLIT           STRING    = Array[3, 2]
;   IDL> print, split
;   2006 Jan 10
;   2009 Feb 10
; note that this is a pretty dumb piece of code and highly susceptible
; to breaking

    nstr = n_elements(str)
    if (nstr eq 0) or (n_elements(key) eq 0) then begin
       doc_library, 'im_strsplit'
       return, -1
    endif
    split1 = strsplit(str[0],key,_extra=extra)
    split = strarr(n_elements(split1),nstr)
    for ii = 0L, nstr-1 do split[*,ii] = strsplit(str[ii],key,_extra=extra)
    
return, split
end
    
