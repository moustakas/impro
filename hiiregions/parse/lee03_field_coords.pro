pro lee03_field_coords
; jm06oct23nyu

; ---------------------------------------------------------------------------
; HoII - coordinates from Hodge et al. 1994
; ---------------------------------------------------------------------------
    
    ra = '08:13:52.45' & dec = '70:52:34.2' ; (B1950/NED)

    hii = ['H1-3','H4','H5','H6-8','H9']
    hiira = ['08:14:14.8','08:14:16.1','08:14:16.6','08:14:16.9','08:14:16.7']
    hiidec = ['70:51:22','70:51:34','70:51:45','70:52:15','70:52:25']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    niceprint, replicate('HoII ',5), hii, raoffset, decoffset
    print
    
return
end
    
