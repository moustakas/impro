function measure_wrapper, name, stars, pscale

    spawn, ['sex J1426p3529_Bw.fits '+name+' -c '+'dwfs.sex'+' -CATALOG_NAME temp.cat'+$
            ' -MAG_ZEROPOINT 32.738 -GAIN 24.5 -SEEING_FWHM 1.21']

; read the catalog

    readsex, 'temp.cat', cathead, catdata

    xpix = sexarr(catdata,cathead,'X_IMAGE')
    ypix = sexarr(catdata,cathead,'Y_IMAGE')
    fwhmobj = sexarr(catdata,cathead,'FWHM_IMAGE')

; read the starlist

    srcor, stars[0,*], stars[1,*], xpix, ypix, 0.5, indx1, indx2, /option ; match starlists
    fwhm = median(fwhmobj[indx2]*pscale)
    
return, fwhm
end

pro pick_sigma
; jm01may28uofa

    path = '/home/ioannis/sirtf/ndwfs/'

    root = 'J1426p3529_'
    filter = 'I'
    pscale = 0.258
    
    imfwhm = 1.210  ; arcsec
    newfwhm = 1.615 ; arcsec (degrade to this FWHM)

    sig = sqrt((newfwhm/2.35/pscale)^2 - (imfwhm/2.35/pscale)^2) ; pixels (nominal smoothing size)
    sigmas = [0.6*sig, 0.8*sig, sig]

    names = root+filter+['_image0.fits','_image1.fits','_image2.fits']
    fwhm = fltarr(3)

    readfast, root+'starlist.txt', stars ; known image stars

    for j = 0L, 2L do begin
    
       dum = measure_wrapper(names[j],stars,pscale)
       fwhm[j] = dum

    endfor

    plotsym, 0, 1, /fill
    window, 0, xs=400, ys=400
    plot, fwhm, sigmas, ps=8, xsty=3, ysty=3, xtit='Resulting FWHM (arcsec)', ytit='Gaussian Sigma (pixels)'

; fit a line and read off the correct sigma

    coef = linfit(sigmas,fwhm)
    goodsig = poly(newfwhm,coef)

    print, 'Final sigma (pixels) for '+filter+' image = '+strn(goodsig)+'.' ; pixels
    
stop

return
end
