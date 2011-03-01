pro view_images
; jm01sep1uofa

    flist = findfile()
    fgood = where(strmatch(flist,'*.fit') eq 1B,fcount) ; FITS files
    flist = flist[fgood]

    fstart = 0L
    fend = 3L
    
    cube = djs_readmany(flist[fstart:fend])

    slaffi, cube, rebin=0.25, imrange=[18000.0,19000.0]

return
end
