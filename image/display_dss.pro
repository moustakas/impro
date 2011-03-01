pro display_dss, image, header=header, write=write, outpath=outpath, $
  title=title, outfile=outfile
; jm02dec11uofa

; IMAGE can either be a FITS file name or a two-dimensional image
; this routine only does one image at a time
    
    nimage = n_elements(image)
    if nimage eq 0L then begin
       print, 'Syntax - display_dss, image'
       return
    endif
    
    imtype = size(image,/type)
    if n_elements(outpath) eq 0L then outpath = cwd()
    if n_elements(outfile) eq 0L then outfile = 'dss_object.png'

    if not keyword_set(write) then window, 0, xs=650, ys=650

    if (imtype eq 7L) then begin

       if file_test(image,/regular) eq 1L then begin

          splog, 'Reading '+image+'.'
          im = readfits(image,header,/silent)

          imsize = size(im,/dimension)
          xsize = imsize[0]
          ysize = imsize[1]

          xmid = float(xsize-1.0)/2.0 & ymid = float(ysize-1.0)/2.0
          pixscale = abs(sxpar(header,'CDELT1'))*3600.0 ; pixel scale [arcsec/pixel]
          
          xaxis = (findgen(xsize)-xmid)*pixscale ;/60.0
          yaxis = (findgen(ysize)-ymid)*pixscale ;/60.0

          xtitle = 'x (arcsec)' & ytitle = 'y (arcsec)'
          
       endif else begin

          splog, 'Image '+image+' not found...returning.'
          return

       endelse

    endif else begin

       im = image

       imsize = size(im,/dimension)
       xsize = imsize[0]
       ysize = imsize[1]

       xaxis = findgen(xsize) & yaxis = findgen(ysize)
       xtitle = 'x (pixel)' & ytitle = 'y (pixel)'
       
    endelse 

    djs_iterstat, im, sigma=imsig, median=immed, sigrej=3.0
;   immin = immed-3*imsig
;   imtop = immed+8*imsig
    
    imdisp, imgscl(im,min=immin,top=imtop,/log), /erase, /axis, /negative, yminor=3, position=[0.15,0.05,0.95,1.0], $
      yrange=minmax(yaxis), xrange=minmax(xaxis), charsize=2.0, charthick=2.0, $
      xtitle=xtitle, ytitle=ytitle, xthick=2.0, ythick=2.0, title=title

; overplot the cardinal directions

    if n_elements(header) ne 0L then $
      arrows, header, 0.8*max(xaxis), 0.7*max(yaxis), /data, thick=4.0
    
; grab the image
    
    if keyword_set(write) then begin
       outim = tvrd()
       splog, 'Writing '+outpath+outfile+'.'
       write_png, outpath+outfile, outim
    endif else begin
;      print, 'Press any key to continue.'
;      cc = strupcase(get_kbrd(1))
    endelse
       
return
end
