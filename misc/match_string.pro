function match_string, starname, stdfiles, index=index, findex=findex, $
  exact=exact, silent=silent
; jm01oct18uofa
; match the standard star data file names to the observed standard
; star name given in the header
; jm02jun25uofa - added exact keyword
; jm02nov12uofa - added FINDEX keyword (first match)    
; jm03may22uofa - fixed case if NAME is ''
; jm06jan19uofa - check for STRMATCH wildcard characters    
    
    nstarname = n_elements(starname)
    if nstarname gt 1L then begin
       findex = lonarr(nstarname)
       for istar = 0L, nstarname-1L do begin
          datfile1 = match_string(starname[istar],stdfiles,index=index1,$
            findex=findex1,exact=exact,silent=silent)
          if istar eq 0L then begin
             index = index1 
             datfile = datfile1 
             findex = findex1
          endif else begin
             index = [index,index1]
             datfile = [datfile,datfile1]
             findex = [findex,findex1]
          endelse
       endfor
       return, datfile
    endif
    
    name = strlowcase(strcompress(starname,/remove))
    name = repstr(repstr(name,'[',''),']')
    len = strlen(name)
    newname = strarr(len>1)
    for i = 0L, len-1L do newname[i] = strmid(name,i,1)

    std = strlowcase(strcompress(stdfiles,/remove))
    std = repstr(repstr(std,'[',''),']')

    if keyword_set(exact) then begin

       index = where(strmatch(std,name) eq 1B,nmatch)
       if (nmatch eq 0L) then begin
          if not keyword_set(silent) then message, 'Exact match not found for '+strupcase(name)+'.', /info
          datfile = ''
          findex = -1L
       endif else begin
          datfile = stdfiles[index]
          findex = index[0]
       endelse

    endif else begin
    
; cycle on each letter in the star name
    
       j = -1L
       repeat begin
          j = j+1L
;         match = where(strcmp(strjoin(newname[0:j]),stdfiles,j+1L,/fold_case) eq 1B,nmatch)
          index = where(strcmp(strjoin(newname[0:j]),std,j+1L,/fold_case) eq 1B,nmatch)
       endrep until (nmatch eq 1B) or (j eq len-1L)

       if (nmatch eq 0B) and (j eq len-1L) then datfile = '' else begin
          datfile = (stdfiles[index])[0]
          findex = index[0]
       endelse
;      message, 'No data file matches!' else datfile = (stdfiles[index])[0]

    endelse

return, datfile
end

