pro pagel79_coords
; jm06dec14nyu - coordinates and numbering from Deharveng et
;                al. 1988 (D88)

    ra = '00:52:31.35' & dec = '-37:57:15.9' ; (B1950; D88)

    phii = ['P1','P5','W21','W23',$
      'P4','P7']
    hii = ['D77','D53A','D137B','D137D',$
      'D159','D6']
    hiira = ['00:52:28.57','00:52:21.67','00:52:51.06','00:52:49.96',$
      '00:53:12.3','00:51:54.49']
    hiidec = ['-37:54:34.9','-37:59:22.8','-37:57:37.9','-37:57:29.4',$
      '-37:59:25.6','-37:51:08.2']

    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'NGC300 ', phii[i], hii[i], $
      hiira[i], hiidec[i], raoffset[i], decoffset[i], format='(A3,A5,A8,2A12,2I7)'
    
Return
end
    
