pro pixon_image

    datapath = '/home/ioannis/sirtf/simulations/'
    imname = 'sim24.fits'
    psfpath = datapath+'psfs/mips_024_newresol_02.00.fits'
    
    image = readfits(datapath+imname,head)
    psf = readfits(psfpath,psfhead)
    
; examine the noise properties of the image by setting pxnmap equal to
; an array of 1's

    imsize = size(image,/dimension)
    xsize = imsize[0]
    ysize = imsize[1]

    pxnmap = fltarr(xsize,ysize)+1.0

    pxn, image, psf, pxnmap, newim, residuals, wsc=0.5
    
    


return
end
