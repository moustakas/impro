;+
; NAME:
;   IM_FIND_NANS
; PURPOSE:
;   Locate NANs in an input data structure.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Oct 11, UCSD
;-

pro im_find_nans, str

    nobj = n_elements(str)

    tags = tag_names(str[0])
    ntag = n_elements(tags)
    for ii = 0, ntag-1 do begin
       len = n_elements(str[0].(ii))
       if (size(str[0].(ii),/type) ne 7) then begin
          inf = where(finite(str.(ii)) eq 0,ninf)
          if (ninf ne 0) then begin
             badobj = inf / len
             badobj = badobj[uniq(badobj,sort(badobj))]
             splog, 'Tag '+tags[ii]+', infinities: '+$
               strjoin(strtrim(badobj,2),',')
          endif
       endif
    endfor
    
return
end
    
