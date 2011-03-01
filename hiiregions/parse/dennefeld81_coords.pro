pro dennefeld81_coords
; jm06dec13nyu - M31 coordinates from Baade & Arp 1964
    
    ra = '00:37:16.94' & dec = '40:43:14.9' ; (B1900; BA64)

    hii = ['BA204','BA289','BA359','BA429','BA479',$
      'BA500','BA666','BA677']
    hiira = ['00:40:13.92','00:36:03.30','00:34:19.52','00:35:34.78','00:33:44.94',$
      '00:32:10.00','00:41:04.54','00:41:01.37']
    hiidec = ['41:20:11.9','40:18:03.2','39:48:00.0','40:31:12.8','40:04:23.5',$
      '39:27:30.2','41:38:35.0','41:40:58.5']
    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'M31 ', hii[i], raoffset[i], decoffset[i]
    
Return
end
    
