pro hunter99_coords
; jm06oct23nyu

; ---------------------------------------------------------------------------
; DDO 53
; ---------------------------------------------------------------------------
    
    ra = '08:34:07.2' & dec = '66:10:54' ; (J2000/NED)

    hii = ['H1','H2','H5','H6']
    hiira = ['08:34:07.2','08:34:08.7','08:34:11.5','08:34:03.8']
    hiidec = ['66:10:54','66:10:52','66:10:11','66:10:37']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, replicate('DDO 53 ',4), hii, raoffset, decoffset
    print
    
return
end
    
