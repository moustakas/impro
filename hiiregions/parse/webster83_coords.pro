pro webster83_coords
; jm06dec14nyu - coordinates and numbering from Deharveng et
;                al. 1988 (D88)

; WS9 and WS12 too faint to be in the D88 catalog; same coordinates
; given for WS14 and WS15
    
    ra = '00:52:31.35' & dec = '-37:57:15.9' ; (B1950; D88)

    wshii = ['WS1','WS2','WS3','WS4','WS5',$
      'WS6','WS7','WS8','WS10','WS11',$
      'WS13','WS14','WS15','WS16']
    hii = ['D100','D76B','D64','D45','D53B',$
      'D119A','D137A','D137C','D134','D24',$
      'D142','D6','D6','D5']
    hiira = ['00:52:35.1','00:52:29.01','00:52:24.8','00:52:18.72','00:52:21.1',$
      '00:52:42.08','00:52:51.15','00:52:52.22','00:52:49.6','00:52:06.95',$
      '00:52:54.02','00:51:54.49','00:51:54.49','00:51:54.29']
    hiidec = ['-37:57:24.2','-37:56:30.3','-37:56:33.9','-37:57:06.9','-37:59:23.7',$
      '-37:59:33.2','-37:57:51.8','-37:57:51.6','-38:00:37.5','-37:57:48.2',$
      '-38:00:33.7','-37:51:08.2','-37:51:08.2','-37:50:48.1']

    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'NGC300 ', wshii[i], hii[i], $
      hiira[i], hiidec[i], raoffset[i], decoffset[i], format='(A3,A5,A8,2A12,2I7)'
    
Return
end
    
