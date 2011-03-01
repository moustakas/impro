function im_today
; jm01oct18uofa
; jm05apr08uofa - remove blank spaces
; return today's date
return, strmid(systime(),20)+' '+strjoin(strsplit(strmid(systime(),4,12),' ',/extract),' ')
end
