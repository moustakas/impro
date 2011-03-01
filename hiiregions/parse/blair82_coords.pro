pro blair82_coords
; jm06dec13nyu - M31 coordinates from Baade & Arp 1964

; coords for BA1/2 taken to be BA1, and for BA487/8 to be BA487     
    
    ra = '00:37:16.94' & dec = '40:43:14.9' ; (B1900; BA64)

    hii = ['BA75','BA423','BA289','BA1/2','BA577',$
      'BA379','BA668','BA676','BA684',$
      'BA487/8','BA381']
    hiira = ['00:38:28.37','00:36:03.52','00:36:03.30','00:28:05.84','00:38:33.08',$
      '00:33:53.72','00:41:04.52','00:41:02.78','00:40:46.98',$
      '00:32:21.03','00:33:12.18']
    hiidec = ['40:53:37.8','40:32:07.0','40:18:03.2','40:37:00.0','41:16:17.8',$
      '39:49:00.0','41:38:38.6','41:39:49.1','41:52:07.7',$
      '39:52:14.1','39:19:54.1']

    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'M31 ', hii[i], hiira[i], hiidec[i], $
      raoffset[i], decoffset[i], format='(A3,A8,2A12,2I7)'
    
Return
end
    
