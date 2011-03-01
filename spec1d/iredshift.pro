;+
; NAME:
;	IREDSHIFT()
;
; PURPOSE:
;       Derive the redshift of a galaxy using PCA eigentemplates and
;       chi-squared fitting.
;
; CALLING SEQUENCE:
;       zans = iredshift(speclist,datapath=,zmin=,zmax=$
;          npoly=,obsname=,eigendir=,eigenfile=,_extra=extra,$
;          /update,/heliocor,/doplot,)
;
; INPUTS:
;	speclist - FITS list of spectra
;
; OPTIONAL INPUTS:
;	datapath  - I/O data path
;	zmin      - minimum redshift to consider (can be negative for 
;                   blushifts) [default -0.01] 
;	zmax      - maximum redshift to consider [default 0.05]
;	npoly     - number of polynomial "eigenspectra" to append to 
;                   the PCA eigentemplates (default 2)
;	obsname   - name of the observatory when HELIOCOR=1 (default 
;                   KPNO)
;       eigenfile - fitting eigentemplates (default to the most recent
;                   "spEigenGal-*") 
;       eigendir  - path name to EIGENFILE (default to
;                   "${IDLSPEC2D_DIR}/templates") 
;	extra     - extra inputs for ZFIND()
;	
; KEYWORD PARAMETERS:
;	update   - updated header corresponding to each object in
;                  SPECLIST with the redshift information
;	heliocr  - modify the derived redshift to a heliocentric frame
;                  of reference using the RA, DEC, EPOCH, and JD
;                  header information
;	doplot   - generate a plot of the spectrum and best-fitting
;                  linear combination of eigentemplates
;
; OUTPUTS:
;	zans     - redshift computation output from ZFIND()
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;	The heliocentric correction to the redshift is placed in the
;	header, but not returned to the user.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	CWD(), SPLOG, DJS_FILEPATH(), READFITS(), SXPAR(),
;	POLY_ARRAY(), RD1DSPEC(), SXADDPAR, COMBINE1FIBER, DJS_MEAN(),
;	ZFIND(), IM_HMS2DEC(), OBSERVATORY, HELIOCENTRIC(), DJS_MODFITS,
;	ICLEANUP, DJS_PLOT, DJS_OPLOT
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 20, U of A
;       jm03jan6uofa - minor updates and modifications
;
; Copyright (C) 2002-2003, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function iredshift, speclist, datapath=datapath, zmin=zmin, zmax=zmax, $
  npoly=npoly, obsname=obsname, eigendir=eigendir, eigenfile=eigenfile, $
  _extra=extra, update=update, heliocor=heliocor, doplot=doplot

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - zans = iredshift(speclist,datapath=,zmin=,zmax=$'
       print, '   npoly=,obsname=,eigendir=,eigenfile=,_extra=extra,$'
       print, '   /update,/heliocor,/doplot,)'
       return, -1
    endif

    light = 2.99792458D5 ; speed of light [km/s]

    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(npoly) eq 0L then npoly = 2L

; redshift search radius

    if n_elements(zmin) eq 0L then zmin = -0.01
    if n_elements(zmax) eq 0L then zmax = +0.05
    
    if n_elements(obsname) eq 0L then obsname = 'KPNO' ; observatory structure

; if multiple spectra have been passed, call this routine recursively

    if nspec gt 1L then begin
       t0 = systime(/seconds)
       for k = 0L, nspec-1L do begin

          zans1 = iredshift(speclist[k],datapath=datapath,zmin=zmin,zmax=zmax,$
            npoly=npoly,obsname=obsname,eigendir=eigendir,eigenfile=eigenfile,$
            _extra=extra,update=update,heliocor=heliocor,doplot=doplot)
          if (k eq 0L) then zans = zans1 else zans = [[zans],[zans1]]
          splog, 'Object #', k, '  Elap time=', systime(/seconds)-t0, $
            ' (sec)  z=', zans1[0].z, ' (pix)'
          
       endfor
       print
       return, zans
    endif

; determine the logarithmic wavelength spacing of the eigentemplates

    if n_elements(eigendir) eq 0L then eigendir = $
      concat_dir(getenv('IDLSPEC2D_DIR'),'templates')
    if n_elements(eigenfile) eq 0L then eigenfile = 'spEigenGal-*'
    allfiles = findfile(djs_filepath(eigenfile,root_dir=eigendir),count=ct)
    if (ct eq 0) then message, 'Unable to find EIGENFILE matching '+eigenfile+'.'
    thisfile = allfiles[(reverse(sort(allfiles)))[0]]

    eigenflux = readfits(thisfile,shdr)
    naxis = sxpar(shdr,'NAXIS1')
    coeff0 = sxpar(shdr,'COEFF0')
    coeff1 = sxpar(shdr,'COEFF1')
    
    if (npoly ne 0L) then eigenflux = [[eigenflux],[poly_array(naxis,npoly)]]
    
    lam = 10D^(coeff0+coeff1*findgen(naxis)) ; logarithmic wavelength vector
    newloglam = alog10(lam)
          
    scube = rd1dspec(speclist,datapath=datapath) ; read in the galaxy spectrum

    objflux = scube.spec        ; object flux [erg/s/cm^2/A]
    objsig = scube.sigspec      ; object sigma
    wave = scube.wave           ; wavelength [A]
    logwave = alog10(wave)

    objivar = objsig*0.0+1.0
    good = where(objsig ne 0.0,ngood)
    if (ngood ne 0L) then objivar[good] = 1.0/objsig[good]^2.0 ; object inverse variance

    header = scube.header
    hdr = header
    sxaddpar, hdr, 'COEFF0', coeff0
    sxaddpar, hdr, 'COEFF1', coeff1
       
; resample the galaxy spectrum into log-lamba, at the spacing of the
; eigenflux log-lambda vector
    
    combine1fiber, logwave, objflux, objivar, newloglam=newloglam, $
      newflux=newflux, newivar=newivar

    fluxnorm = double(djs_mean(newflux)) ; normalization constants
    ivarnorm = double(djs_mean(newivar))

    if keyword_set(doplot) and (!d.window ne 0L) then window, 0L, xs=450, ys=450
    zans = zfind(newflux/fluxnorm,newivar/ivarnorm,hdr=hdr,npoly=0,$
      starflux=eigenflux,starloglam0=coeff0,stardloglam0=coeff1,$
      zmin=zmin,zmax=zmax,doplot=doplot,_extra=extra)
      
; compute the heliocentric correction

    if keyword_set(heliocor) then begin

       hmsra = sxpar(header,'RA',count=racount)
       dmsdec = sxpar(header,'DEC',count=deccount)
       epoch = sxpar(header,'EPOCH',count=ecount)
       jd = sxpar(header,'JD',count=jdcount)

       if (racount+deccount+ecount+jdcount)[0] ne 4L then begin

          splog, 'HELIOCOR: Insufficient header information in '+speclist+'.'
          z_helio = zans[0].z

       endif else begin
       
          observatory, obsname, obs 

          ra = 15D*im_hms2dec(hmsra)
          dec = im_hms2dec(dmsdec)
       
          vcorr = heliocentric(ra,dec,epoch,jd=jd,longitude=(360-obs.longitude),$
            latitude=obs.latitude,altitude=obs.altitude)

          z_helio = zans[0].z + vcorr[0]/light

       endelse
          
    endif else z_helio = zans[0].z

; compute the velocity dispersion

;   vdans = vdispfit(newflux/fluxnorm,newivar/ivarnorm,hdr=hdr,newloglam,$
;     npoly=npoly,zobj=zans[0,*].z,eigendir=eigendir,eigenfile='spEigenElodie.fits',$
;     columns=lindgen(24),yfit=dispflux,/doplot)

; update the header

    if keyword_set(update) then begin

       splog, 'Updating the header for '+datapath+speclist+'.'
       newhead = header
       sxaddpar, newhead, 'Z', float(z_helio), ' redshift', before='HISTORY'
       sxaddpar, newhead, 'Z_ERR', zans.z_err, ' redshift error', before='HISTORY'

       djs_modfits, datapath+speclist, 0, newhead

    endif

; plot the spectrum and the best fit.  this plot doesn't work if the
; COLUMNS keyword is used

    if keyword_set(doplot) and (not keyword_set(columns)) then begin

       if (!d.window ne 2L) then window, 2, xs=450, ys=450

       if n_elements(columns) ne 0L then neigen = n_elements(columns) else $
         neigen = (size(eigenflux,/dimension))[1]
       bestspec = newflux*0.0
       for j = 0L, neigen-1L do bestspec = bestspec + zans.theta[j]*eigenflux[*,j]
       
       djs_plot, lam/(1+z_helio), newflux/fluxnorm, ps=10, $
         xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
         ythick=2.0, xtitle='Wavelength', ytitle='Relative Flux', $
         color='red'
       djs_oplot, lam, bestspec, ps=10, line=0, color='green'
;      cc = get_kbrd(1)
       
    endif

return, zans
end
