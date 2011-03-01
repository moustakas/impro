pro miller96_coords
; jm06oct23nyu

; ---------------------------------------------------------------------------
; Holmberg I
; ---------------------------------------------------------------------------
    
    ra = '09:36:04.2' & dec = '71:24:33' ; (B1950/NED)

    hii = ['MH25']
    hiira = ['09:36:10.0']
    hiidec = ['71:23:51']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, 'Holmberg I', hii, raoffset, decoffset
    print
    
; ---------------------------------------------------------------------------
; M81DwB
; ---------------------------------------------------------------------------
    
    ra = '10:01:25.09' & dec = '70:36:27.2' ; (B1950/NED)

    hii = ['MH1']
    hiira = ['10:01:23.2']
    hiidec = ['70:36:40']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, 'M81 DwB ', hii, raoffset, decoffset
    print
    
; ---------------------------------------------------------------------------
; IC 2574
; ---------------------------------------------------------------------------
    
    ra = '10:24:42.35' & dec = '68:40:03.5' ; (B1950/NED)

    hii = ['MH167','MH198','MH242']
    hiira = ['10:25:02.5','10:25:08.6','10:25:18.0']
    hiidec = ['68:43:47','68:43:48','68:43:50']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, replicate('IC2574 ',3), hii, raoffset, decoffset
    print

; test:    
    
;   print
;   cosd = cos(im_hms2dec('68:24:43.7')*!dtor)
;   rr = im_hms2dec('10:28:23.48')*15.0+194.52D/3600.0D/cosd
;   print, im_dec2hms(rr/15.0D)
    
Return
end
    
