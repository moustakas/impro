; ----------------------------------------------------------------------

function str2time,str_time

str_time = strtrim(str_time,2)
if (strmid(str_time,0,1) EQ "'") then begin
   str_time = strmid(str_time,1,strlen(str_time)-1)
endif
time = float(strmid(str_time,0,2))
time = time + float(strmid(str_time,3,2))/60.0
time = time + float(strmid(str_time,6,2))/3600.0

return,time

end

; ----------------------------------------------------------------------

function time2str,time

hours = fix(time)
min = round(60.0*(time - hours))
if (min LE 9) then begin
    add_zero = '0'
endif else begin
    add_zero = ""
endelse

str_time = strtrim(string(hours),2) + ":" + add_zero + strtrim(string(min),2)

return,str_time

end

; ----------------------------------------------------------------------
