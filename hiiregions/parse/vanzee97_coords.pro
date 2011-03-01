pro vanzee97_coords
; jm06oct23nyu

; ---------------------------------------------------------------------------
; DDO 154 = UGC 08024
; ---------------------------------------------------------------------------
    
    ra = '12:51:39.37' & dec = '27:25:14.2' ; (B1950/NED)

    hii = ['H1','H2']
    hiira = ['12:51:44.3','12:51:41.6']
    hiidec = ['27:25:26','27:25:26']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, replicate('DDO 154 ',2), hii, raoffset, decoffset
    print
    
return
end
    
