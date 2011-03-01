pro dodorico83_coords
; jm06dec14nyu - coordinates and numbering from Deharveng et
;                al. 1988 (D88)

; same coordinates used for D137D (W21 & W23)
    
    ra = '00:52:31.35' & dec = '-37:57:15.9' ; (B1950; D88)

    dhii = ['H1','H2','H3','H4','H5','H6',$
      'H7','H8','H9','H10','H11',$
      'H12','H13','H14','H15']
    hii = ['D6','D24','D45','D52','D53C','D53A',$
      'D77','D79','D109','D115','D118A',$
      'D119A','D137D','D137A','D137C']
    hiira = ['00:51:54.49','00:52:06.95','00:52:18.72','00:52:21.13','00:52:21.08','00:52:21.67',$
      '00:52:28.57','00:52:29.41','00:52:38.68','00:52:40.93','00:52:41.97',$
      '00:52:42.08','00:52:49.96','00:52:51.15','00:52:52.22']
    hiidec = ['-37:51:08.2','-37:57:48.2','-37:57:06.9','-37:56:16.8','-37:59:10.4','-37:59:22.8',$
      '-37:54:34.9','-37:54:37.9','-37:56:48.0','-37:54:41.9','-37:59:01.2',$
      '-37:59:33.2','-37:57:29.4','-37:57:51.8','-37:57:51.6']

    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'NGC300 ', dhii[i], hii[i], $
      hiira[i], hiidec[i], raoffset[i], decoffset[i], format='(A3,A5,A8,2A12,2I7)'
    
Return
end
    
