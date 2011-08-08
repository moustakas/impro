;+
; NAME:
;   FITSDATE2JD()
; PURPOSE:
;   Convert FITS format date to julian date.
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 Apr 23, U of A
;-

function fitsdate2jd, date

    ndate = n_elements(date)
    if ndate eq 0L then begin
       print, 'Syntax - date = fitsdate2jd(date)'
       return, -1
    endif

    jd = dblarr(ndate)
    for i = 0L, ndate-1L do begin
    
       d = strsplit(date[i],'-',/extract)
       jd[i] = julday(d[1],d[2],d[0])

    endfor

    if ndate eq 1L then jd = jd[0]

return, jd
end
