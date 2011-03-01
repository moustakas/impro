pro kennicutt01_coords
; jm06oct23nyu

; ---------------------------------------------------------------------------
; DDO 154
; ---------------------------------------------------------------------------
    
    ra = '12:54:05.2' & dec = '27:08:59' ; (J2000/NED)

    hii = ['H2','H3']
    hiira = ['12:54:03.9','12:54:03.9']
    hiidec = ['27:09:05','27:08:26']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, replicate('DDO 154 ',2), hii, raoffset, decoffset
    print
    
return
end
    
