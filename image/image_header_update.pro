pro image_header_update

    path = pisces_path()
    pushd, path

    flist = findfile('*.fits',count=fcount)

    obsname = 'MMT'
    obsinfo = {observatory: obsname, name:'Multiple Mirror Telescope', $
               longitude: ten([110,53,4.4D]), latitude: ten([31,41,19.6D]), $
               altitude: 2606D, tz: 7.0}
 
    ra = '14:26:00.8'
    dec = '+35:31:31'
    epoch = '2000.0'
    dateobs = '2001-03-14' ; date of observation at start of night
    ymd = strsplit(dateobs,'-',/extract)
    
    for j = 0L, fcount-1L do begin

;      h = headfits(flist[j])
       data = readfits(flist[j],h,/silent)
       h = h[where(strcompress(h,/remove) ne '')]
       
       time = sxpar(h,'TIME')
       hms = float(strsplit(time,':',/extract))
       hms[0] = hms[0] + 12.0
       
       if long(hms[0])/24L then begin
          hms[0] = hms[0] - 24.0
          ymd[2] = ymd[2] + 1   ; advance the date by one day
       endif
       
       ymd = strcompress(ymd,/remove)
       newdate = strjoin(ymd,'-')

       dtime = ten(hms) ; decimal time
       newtime = strjoin(strsplit(im_dec2hms(dtime),' ',/extract),':')
       
; compute the local sidereal time
       
       ct2lst, dlst, (360.0-obsinfo.longitude), obsinfo.tz, dtime, ymd[2], ymd[1], ymd[0]
       lst = strjoin(strsplit(im_dec2hms(dlst),' ',/extract),':')
       
       ut = im_dec2hms((dtime+obsinfo.tz) mod 24.0)              ; universal time
       newut = strjoin(strsplit(ut,' ',/extract),':')
       jd = julday(ymd[1],ymd[2],ymd[0],hms[0],hms[1],hms[2]) ; julian date

       lha = dlst*15.0D - 15.0*im_hms2dec(ra) ; local hour angle
       
; compute the airmass, the zenith distance, and the parallactic angle

       airmass = compute_airmass(obsinfo.latitude,im_hms2dec(dec),lha,zd=zd)
       parangle = compute_parangle(obsinfo.latitude,im_hms2dec(dec),lha)
       
; update the header

       sxdelpar, h, 'TIME'
       sxdelpar, h, 'LHA'
       sxdelpar, h, 'SECZ'

       sxaddpar, h, 'RA', ra, 'Right ascension [HMS]'
       sxaddpar, h, 'DEC', dec, 'Declination [DMS]'
       sxaddpar, h, 'EPOCH', epoch, after='DEC'
       sxaddpar, h, 'MST', newtime, 'Mountain standard time', after='EPOCH'
       sxaddpar, h, 'DATE-OBS', newdate, 'Mountain standard date', after='MST'
       sxaddpar, h, 'UT', newut, 'Universal time', after='MST'
       sxaddpar, h, 'JD', jd, 'Julian date', after='UT', format='(F22.10)'
       sxaddpar, h, 'LST', lst, 'lst', after='JD', 'Local sidereal time'
       sxaddpar, h, 'HA', lha, after='LST', 'Hour angle [degrees]', format='(F6.2)'
       sxaddpar, h, 'AIRMASS', airmass, after='HA', 'Airmass', format='(F6.4)'
       sxaddpar, h, 'ZD', zd, after='AIRMASS', 'Zenith distance [degrees]', format='(F6.2)'
       sxaddpar, h, 'PA', parangle, after='ZD', 'Parallactic angle [degrees]', format='(F6.2)'

; write out the new headers

;      djs_modfits, flist[j], data, h
;      mwrfits, data, flist[j], h
       writefits, flist[j], data, h

    endfor

    popd
    
return
end    
