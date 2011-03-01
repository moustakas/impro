function jd2fitsdate, jd, noyear=noyear
; jm02oct17uofa
; convert julian date to a FITS format date    

    njd = n_elements(jd)
    if njd eq 0L then begin
       print, 'Syntax - date = jd2fitsdate(jd)'
       return, -1
    endif

    daycnv, jd, yr, mn, day, hr

    yr = string(yr,format='(I4)')
    mn = string(mn,format='(I2.2)')
    day = string(day,format='(I2.2)')

    date = strarr(njd)
    if keyword_set(noyear) then $
      for i = 0L, njd-1L do date[i] = strjoin([mn[i],day[i]],'-') else $
      for i = 0L, njd-1L do date[i] = strjoin([yr[i],mn[i],day[i]],'-')

    if njd eq 1L then date = date[0]

return, date
end
