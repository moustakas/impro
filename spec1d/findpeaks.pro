pro findpeaks, spec, x, y, width
; jm01jul20uofa
; peak-finding algorithm.  basically, loop on each pixel and demand
; that a peak is defined as having all pixel values on either side of
; the peak within a "box" radius to be below the peak value.  we also
; require those pixel values to be monotonically increasing or
; decreasing.  

; add a gaussian keyword or a flux-weighted keyword or whatever

; need to use the variance!
    
; read in an emission (arclamp) spectrum and the variance spectrum

; need some strong criteria for acceptance or rejection of a line,
; besides just a FWHM check
    
;   spec = readfits('/home/ioannis/kennicutt/arc_spec.fits',head,/silent)
;   vspec = readfits('/home/ioannis/kennicutt/arc_vmap.fits',/silent)

    box = 3L ; boxsize radius (pixels)
    npix = size(spec,/n_elements)
    specaxis = findgen(npix)
    
;   window, 0, xs=450, ys=450
;   plot, specaxis, spec, color=5, xsty=3
;   window, 2, xs=450, ys=450
    plotsym, 0, 2.0, thick=2.5
    colortable2, c

; define a continuum level to define "strong" emission lines

;   clipsig = 2.0
;   djs_iterstat, spec, sigrej=clipsig, mean=smean, sigma=ssigma, mask=specmask

    xpeaks = ptr_new()
    ypeaks = ptr_new()
    fwhm = ptr_new()
    
;   user = ''
    for j = box+1L, npix-box-2L do begin

       area = [spec[j-box:j-1L],spec[j+1L:j+box]] ; data area around the current pixel

       below = where(area lt spec[j],nbelow)
;      if ((nbelow eq 2L*box) and specmask[j] and monotonic(area[0L:box-1L]) and monotonic(area[box:2*box-1L])) then begin
       if ((nbelow eq 2L*box) and monotonic(area[0L:box-1L]) and monotonic(area[box:2*box-1L])) then begin

; check the FWHM of the line using interpolation (don't assume the
; line is gaussian)

          halfmax = spec[j]/2.0
          lhs = interpol(specaxis[j-box:j],spec[j-box:j],halfmax)
          rhs = interpol(specaxis[j:j+box],spec[j:j+box],halfmax)
          width = rhs-lhs          

; do a flux-weighted centering

          onedcenter, specaxis[j-box:j+box], spec[j-box:j+box], xcen, ycen
          
          if ptr_valid(fwhm) then fwhm = ptr_new([*fwhm,width]) else fwhm = ptr_new(width)

;         if ptr_valid(xpeaks) then xpeaks = ptr_new([*xpeaks,specaxis[j]]) else xpeaks = ptr_new(specaxis[j])
          if ptr_valid(xpeaks) then xpeaks = ptr_new([*xpeaks,xcen]) else xpeaks = ptr_new(xcen)
;         if ptr_valid(ypeaks) then ypeaks = ptr_new([*ypeaks,spec[j]]) else ypeaks = ptr_new(spec[j])
          if ptr_valid(ypeaks) then ypeaks = ptr_new([*ypeaks,ycen]) else ypeaks = ptr_new(ycen)
         
;         wset, 0
;         plots, [specaxis[j],specaxis[j]], 1.1*[spec[j],spec[j]], ps=8, color=6, /data
;         wset, 2
;          plot, specaxis[(j-2*box)>0L:(j+2*box)<(npix-1L)], spec[(j-2*box)>0L:(j+2*box)<(npix-1L)], color=4, xsty=3
;          oplot, [!x.crange[0],!x.crange[1]], [smean,smean], line=2, color=8
;          oplot, [specaxis[j-box],specaxis[j-box]], [!y.crange[0],!y.crange[1]], line=2, color=7
;          oplot, [specaxis[j+box],specaxis[j+box]], [!y.crange[0],!y.crange[1]], line=2, color=7
;          plots, [specaxis[j],specaxis[j]], 1.1*[spec[j],spec[j]], ps=8, color=6, /data
;
;          splog, j, width, /noname
;          read, user

       endif

    endfor

    width = *fwhm
    x = *xpeaks
    y = *ypeaks

; eliminate bad lines
    
    djs_iterstat, width, sigrej=3.0, mask=fwhmask
    glines = where(fwhmask eq 0B,nglines)
    if nglines eq 0L then message, 'All lines were eliminated based on their FWHM!' else begin
       x = x[glines]
       y = y[glines]
       width = width[glines]
    endelse

; clean up

    if ptr_valid(fwhm) then ptr_free, fwhm
    if ptr_valid(xpeaks) then ptr_free, xpeaks
    if ptr_valid(ypeaks) then ptr_free, ypeaks

return
end

