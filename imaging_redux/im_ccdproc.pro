;+
; NAME:
;   IM_CCDPROC
;
; PURPOSE:
;   Apply standard imaging reductions to a list of images.
;
; INPUTS: 
;   imagefile - FITS file names to reduce
;
; OPTIONAL INPUTS: 
;   outpath    - output path name (default './')
;   prefix     - prefix to append to the output images (default 'r')
;   gain       - gain [electron/ADU]
;   rdnoise    - read noise [electron]
;   biasfile   - FITS file name of the master bias frame
;   darkfile   - FITS file name of the master dark frame
;   flatfile   - FITS file name of the master flat-field
;   badpixfile - FITS file name of the bad pixel mask (can be
;                a scalar, or an array, one for each image)
;   biassec    - overscan region
;   trimsec    - trim region
;   satvalue   - saturation value (default 50000 ADU)
;
; KEYWORD PARAMETERS: 
;   crrej - identify cosmic rays; add them to the bad pixel mask 
;   wfits - write out the reduced FITS file
;
; OUTPUTS: 
;   Individual images, inverse variance (weight) maps, RMS maps,
;   and FLAG maps are written out.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 Nov, NYU - written
;   jm07jannyu - developed more
;   jm07jun20nyu - added CRREJ keyword
;   jm07aug24nyu - added SATVALUE optional input 
;   jm08jul30nyu - more generalized treatment of masked pixels;
;     WEIGHT_SUFFIX optional input removed
;
; Copyright (C) 2006-2008, John Moustakas
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

pro im_ccdproc, imagefile, outpath=outpath, prefix=prefix, gain=gain, rdnoise=rdnoise, $
  biasfile=biasfile, darkfile=darkfile, flatfile=flatfile, badpixfile=badpixfile1, $
  biassec=biassec, trimsec=trimsec, satvalue=satvalue, crrej=crrej, wfits=wfits, $
  grow_satmask=grow_satmask

    nimage = n_elements(imagefile) 
    if (nimage eq 0) then begin
       doc_library, 'im_ccdproc'
       return
    endif

    if (n_elements(outpath) eq 0) then outpath = './'
    if (n_elements(prefix) eq 0) then prefix = 'r'
    if (n_elements(gain) eq 0) then gain = 1.0
    if (n_elements(rdnoise) eq 0) then rdnoise = 0.0
    if (n_elements(satvalue) eq 0) then satvalue = 50000.0 ; [ADU]

    outfile = outpath+prefix+file_basename(imagefile)

    nbadpixfile = n_elements(badpixfile1)
    if (nbadpixfile ne 0) then begin
       if (nbadpixfile eq 1L) then badpixfile = replicate(badpixfile1,nimage) else begin
          if (nbadpixfile ne nimage) then begin
             splog, 'IMAGEFILE and BADPIXFILE dimensions do not match!'
             return
          endif else badpixfile = badpixfile1
       endelse
    endif
    
    t0 = systime(1)
    for ii = 0, nimage-1L do begin
       if (file_test(imagefile[ii]) eq 0) then begin
          splog, 'Image '+imagefile[ii]+' not found'
          return
       endif
       
       splog, 'Processing image '+imagefile[ii]
;      print, format='("Processing image ",I0,"/",I0,a1,$ )', $
;        ii+1, nimage, string(13b)

; read the image, overscan-subtract, and trim
       im = float(readfits(imagefile[ii],hdr,/silent))
       im_overscan, im, biassec=biassec, hdr=hdr, /median, silent=silent
       im_trim, im, trimsec=trimsec, hdr=hdr, silent=silent

; subtract the master bias frame
       proceed = 1
       if (n_elements(biasfile) eq 0) then proceed = 0 else $
         if (file_test(biasfile) eq 0) then proceed = 0

       if (proceed eq 0) then begin
          splog, 'BIASFILE not defined or does not exist.'
          return
       endif
             
       biasframe = mrdfits(biasfile,0,biashdr,/silent)

       headcheck1 = sxpar(biashdr,'OVERSCAN',count=headcount1)
       headcheck2 = sxpar(biashdr,'TRIM',count=headcount2)
       if (headcount1 eq 0) or (headcount2 eq 0) then begin
          splog, 'Master bias frame needs to have been overscan subtracted and trimmed!' 
          return
       endif

       im_biasc, im, biasframe=biasframe, hdr=hdr, silent=silent, $
         biasfile=file_basename(biasfile)
          
; subtract the master dark frame
;      proceed = 1
;      if (n_elements(darkfile) eq 0) then proceed = 0 else $
;    if (file_test(darkfile) eq 0) then proceed = 0
;
;      if (proceed eq 0) then begin
;     splog, 'DARKFILE not defined or does not exist.'
;     return
;      endif
;        
;      darkframe = mrdfits(darkfile,0,darkhdr,/silent)
;
;      headcheck1 = sxpar(darkhdr,'OVERSCAN',count=headcount1)
;      headcheck2 = sxpar(darkhdr,'TRIM',count=headcount2)
;      if (headcount1 eq 0) or (headcount2 eq 0) then begin
;     splog, 'Master dark frame needs to have been overscan subtracted and trimmed!' 
;     return
;      endif
;
;      im_darkc, im, exptime, darkframe=darkframe, darkfile=file_basename(darkfile), $
;    hdr=hdr, silent=silent

; divide by the flat field; make sure the flat-field has been trimmed
; to the same dimensions as the data image
       proceed = 1
       if (n_elements(flatfile) eq 0) then proceed = 0 else $
         if (file_test(flatfile) eq 0) then proceed = 0

       if (proceed eq 0) then begin
          splog, 'FLATFILE not defined or does not exist.' 
          return
       endif
       
       flatfield = mrdfits(flatfile,0,flathdr,/silent)

       headcheck = sxpar(flathdr,'TRIM',count=headcount)
       if (headcount eq 0) then begin
          splog, 'Master flat-field needs to have been trimmed!' 
          return
       endif

       varmap = sqrt(im^2.0)/gain + (rdnoise/gain)^2.0 ; variance map [ADU^2]

       im_flatfield, im, varmap=varmap, flatfield=flatfield, $
         silent=silent, flatfile=file_basename(flatfile), hdr=hdr

       invmap = 1.0/(varmap + (varmap eq 0.0))*(varmap ne 0.0) ; [1/ADU^2]
       rmsmap = sqrt(varmap) ; [ADU]
       
; read and/or generate the bad pixel mask
       if (nbadpixfile eq 0) then begin
          morebadpix = long(im*0.0)
       endif else begin
          if (file_test(badpixfile[ii]) eq 0) then begin
             splog, 'BADPIXFILE '+badpixfile[ii]+' not found'
             return
          endif else begin
             splog, 'Using the bad pixel mask, '+file_basename(badpixfile[ii])
             morebadpix = mrdfits(badpixfile[ii],0,badpixhdr,/silent)
             headcheck = sxpar(badpixhdr,'TRIM',count=headcount)
             if (headcount eq 0) then begin
                splog, 'Bad pixel mask '+badpixfile[ii]+' needs to have been trimmed!' 
                return
             endif
          endelse
       endelse

       splog, 'Building the flag map'
       flag = long(im*0.0) or morebadpix

; build a saturated pixel mask; grow the hell out of it, if desired 
       satmask = im ge satvalue
       if keyword_set(grow_satmask) then begin
          satmask = smooth(float(satmask),11,/nan,/edge) gt 0.0
          satmask = smooth(float(satmask),11,/nan,/edge) gt 0.0
       endif

       sat = where(satmask gt 0.0,nsat)
       if (nsat ne 0) then flag[sat] = flag[sat] or $
         im_flagval('PIXMASK','SATURATED')

; identify low-response pixels       
       badflat = where((flatfield le 0.75),nbadflat)
       if (nbadflat ne 0) then flag[badflat] = flag[badflat] or $
         im_flagval('PIXMASK','FLATFIELD_0.75')

; reset the "no data" bad pixel values to just "no data", since many
; of the pixels that may have been flagged as FLATFIELD_0.75 are
; likely in those "no data" regions
       nodata = where((flag and im_flagval('PIXMASK','NODATA') ne 0),nnodata)
       if (nnodata ne 0) then $
         flag[nodata] = im_flagval('PIXMASK','NODATA')
       
;      imbadpix = (im le 0.0) or (im ge satvalue) ; 1=bad, 0=good
;      flatbadpix = (flatfield le 0.05)           ; 1=bad, 0=good
;      badpix = (imbadpix + flatbadpix + morebadpix) eq 0B ; switch! 1=good, 0=bad

;      invmap = 1.0 / (im/gain + (rdnoise/gain)^2.0)
;      invmap = invmap * badpix
;      invmap = invmap < smooth(invmap,11,/nan,/edge)
;      im = im * (invmap gt 0.0)

; reject cosmic rays; first do a simple sky subtraction; also set the
; inverse variance of previously flagged pixels to zero, to prevent
; them from being identified as cosmic rays
       if keyword_set(crrej) then begin
;         splog, 'Rejecting cosmic rays'
          sky, im[where(im gt 0.0)], skymode, skysig, $
            readnoise=rdnoise, /silent
          imnosky = im - skymode
          crinvmap = invmap & maskthese = where(flag gt 0,nmaskthese)
          if (nmaskthese ne 0) then crinvmap[maskthese] = 0
          reject_cr, imnosky, crinvmap, [0.496,0.246], rejects, $
            nrejects=nrejects, c2fudge=c2fudge, niter=10
          splog, 'Identified '+string(nrejects,format='(I0)')+' cosmic rays'
          if (nrejects gt 0) then begin
             flag[rejects] = flag[rejects] or im_flagval('PIXMASK','CR_REJECT')
; could interpolate here
;         invmap[rejects] = 0.0
;         im = im * (invmap gt 0.0)
;         imnosky = imnosky*(invmap gt 0.0)
;         imnosky = djs_maskinterp(imnosky,(invmap le 0),iaxis=0,/const)
          endif
       endif

; now deal with bad pixels       
       
; for cosmetic purposes, set the "no data" regions in the image equal
; to zero; don't mask out saturated pixels, because we want to
; track these with our SE catalogs later
       nodata = where((flag and im_flagval('PIXMASK','NODATA') ne 0),nnodata)
       if (nnodata ne 0) then begin
          im[nodata] = 0.0
          invmap[nodata] = 0.0
          rmsmap[nodata] = 1E30
       endif

; for aesthetic reasons, interpolate over bad columns, cosmic ray
; pixels, and satellite trails; set the inverse variance in these
; pixels to zero and the rms value to a large value
       interp = where((flag and im_flagval('PIXMASK','BAD_COLUMN')) or $
         (flag and im_flagval('PIXMASK','BAD_PIXEL')) or $
;        (flag and im_flagval('PIXMASK','FLATFIELD_0.75')) or $
         (flag and im_flagval('PIXMASK','CR_REJECT')) or $
         (flag and im_flagval('PIXMASK','SATTRAIL')),ninterp)

       if (ninterp ne 0) then begin
          splog, 'Interpolating over bad columns, cosmic rays, etc.'
          mask = byte(im*0.0) & mask[interp] = 1B
          im = djs_maskinterp(im,mask,iaxis=0,/const)
          invmap[interp] = 0.0
          rmsmap[interp] = 1E30
       endif

; write out: (1) image; (2) inverse variance (weight) map; (3) RMS map
; (for SExtractor); and (4) flag map (again for SExtractor)
       if keyword_set(wfits) then begin
          sxdelpar, hdr, 'BSCALE'
          sxdelpar, hdr, 'BZERO'
; image
          splog, 'Writing '+outfile[ii]
          mwrfits, float(im), outfile[ii], hdr, /create
; inverse variance (weight) map
          weightfile = repstr(outfile[ii],'.fits','.weight.fits')
          splog, 'Writing '+weightfile
          mwrfits, float(invmap), weightfile, hdr, /create
; RMS map
          rmsfile = repstr(outfile[ii],'.fits','.rms.fits')
          splog, 'Writing '+rmsfile
          mwrfits, float(rmsmap), rmsfile, hdr, /create
; flag map
          flagfile = repstr(outfile[ii],'.fits','.flag.fits')
          splog, 'Writing '+flagfile
          mwrfits, long(flag), flagfile, hdr, /create
       endif 

    endfor
    splog, 'Total time to process all images ', (systime(1)-t0), ' (seconds)'
       
return    
end


   
