;+
; NAME:
;   IM_TODAY()
; PURPOSE:
;   Return today's date
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Oct 08, U of A
;   jm05apr08uofa - remove blank spaces
;-

function im_today
return, strmid(systime(),20)+' '+strjoin(strsplit(strmid(systime(),4,12),' ',/extract),' ')
end
