;+
; NAME:
;   IM_TODAY()
; PURPOSE:
;   Return today's date
; KEYWORD PARAMETERS:
;   justdate - just return the date, not the time, too
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Oct 08, U of A
;   jm05apr08uofa - remove blank spaces
;-

function im_today, justdate=justdate
    if keyword_set(justdate) then x2 = 6 else x2 = 12
    today = strmid(systime(),20)+' '+strjoin(strsplit(strmid(systime(),4,x2),' ',/extract),' ')
return, today
end
