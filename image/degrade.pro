function padim, image, npix, xsize, ysize, restore=restore
; add 50 pixels of padding to an image before convolution

    info = sstretch(image,/silent)

    bigx = xsize+npix*2L
    bigy = ysize+npix*2L
       
    if keyword_set(restore) then newim = image[npix:bigx-npix-1L,npix:bigy-npix-1L] else begin

       skyb = info.median       ; median sky level

       newim = fltarr(bigx,bigy)
       newim[npix:bigx-npix-1L,npix:bigy-npix-1L] = image
       
       newim[0:npix-1L,*] = skyb
       newim[*,0:npix-1L] = skyb
       newim[bigx-npix-1L:bigx-1L,*] = skyb
       newim[*,bigy-npix-1L:bigy-1L] = skyb
       
    endelse

;   window, 0, xs=450, ys=450
;   display, bytscl(newim,min=info.min,max=info.max)

return, newim
end

pro degrade, image, header, sigma, kernal, fitsname
; jm01may25uofa

; sigma - Gaussian sigma of the convolution kernal
; degrade each of the images to the worst seeing

    xsize = sxpar(header,'NAXIS1')
    Ysize = sxpar(header,'NAXIS2')

    temp = padim(image,50L,xsize,ysize)
    newim = gaussconv(temp,sigma,kernal)

    newimage = padim(newim,50L,xsize,ysize,/restore)

    print, 'Writing '+fitsname+'.'
    writefits, fitsname, newimage, header

return
end
