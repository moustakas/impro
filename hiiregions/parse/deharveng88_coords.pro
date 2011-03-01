pro deharveng88_coords
; jm06dec14nyu - coordinates only needed for D84 and D3; the remaining
;                coordinates were taken from Webster & Smith 1983,
;                d'Odorico et al. 1983, or Pagel et al. 1979

    ra = '00:52:31.35' & dec = '-37:57:15.9' ; (B1950; D88)

    hii = ['D84','D3']
    hiira = ['00:52:30.06','00:51:49.19']
    hiidec = ['-37:55:50.9','-37:53:57.1']

    cosd = cos(im_hms2dec(dec)*!dtor)
    decoffset = (im_hms2dec(hiidec)-im_hms2dec(dec))*3600.0
    raoffset = (im_hms2dec(hiira)-im_hms2dec(ra))*15.0*cosd*3600.0
    for i = 0L, n_elements(hii)-1L do print, 'NGC300 ', hii[i], $
      hiira[i], hiidec[i], raoffset[i], decoffset[i], format='(A3,A8,2A12,2I7)'
    
Return
end
    
