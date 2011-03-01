function scaleimage, image, log=log
; jm02oct3uofa
; mostly taken from ATV
    
;   stats = sstretch(image,npix=1000,/silent)
;   im = (((image-stats.median) < (10*stats.sigma)) > (3*stats.sigma))

    if keyword_set(log) then begin

       maximage = max(image)
       minimage = min(image)

       offset = minimage - (maximage-minimage)*0.01
       immin = alog10(minimage - offset)
       immax = alog10(maximage - offset)

       image = alog10(image-offset)
       
    endif else begin
    
       med = median(image)
       sig = stddev(image)
       immax = (med + (10*sig)) < max(image)
       immin = (med - (2*sig)) > min(image)

    endelse

    image = bytscl(image,min=immin,max=immax,/nan,top=!d.table_size-1)

return, image
end
    
