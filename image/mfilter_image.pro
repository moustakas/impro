;+
; NAME:
;	MFILTER_IMAGE
;
; PURPOSE:
;	Median-filter a simulated SIRTF image.
;
; CALLING SEQUENCE:
;	mfilter_image, fitsname, [box=], [outname=]
;
; INPUTS:
;	fitsname - name of the FITS image
;
; OPTIONAL INPUTS:
;	nbeams   - number of beams:  median filter box size (default 3) 
;	outname  - output name of the filtered image
;	datapath - data path for I/O
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;	READFITS(), SXPAR(), STRN(), WRITEFITS, SXADDHIST
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August 31, U of A
;
; Copyright (C) 2001, John Moustakas
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

pro mfilter_image, fitsname, nbeams=nbeams, outname=outname, datapath=datapath

    if n_params() ne 1L then begin
       print, 'Syntax - mfilter_image, fitsname, [nbeams=], [outname=]'
       return
    endif

    if not keyword_set(datapath) then datapath = './'
    if keyword_set(nbeams) then nbeams = long(nbeams) else nbeams = 3L
    
    if (size(fitsname[0],/type) ne 7L) or (strmatch(fitsname[0],'*.fits') eq 0B) then $
      message, 'Image names must be type string FITS files.'
    nimage = n_elements(fitsname)

    if keyword_set(outname) then begin
       if (size(outname[0],/type) ne 7L) or (strmatch(outname[0],'*.fits') eq 0B) then $
         message, 'Output image names must be type string FITS files.'
       if nimage ne n_elements(outname) then message, 'FITSNAME and OUTNAME must have the same number of elements!'
    endif else begin
       outname = strarr(nimage)
       for k = 0L, nimage-1L do outname[k] = strmid(fitsname[k],0,strpos(fitsname[k],'.fits'))+'_mfilter'+$
         strn(nbeams)+'.fits'
    endelse

    for i = 0L, nimage-1L do begin
       
       print, 'Reading '+strn(fitsname[i])
       image = readfits(datapath+fitsname[i],header,/silent) 

; calculate the beam size
       
       lambda0 = sxpar(header,'WAVE')          ; central wavelength (micron)
       pscale = 3600.0 *sxpar(header,'CDELT1') ; pixel scale (arcsec/pixel)

       beam = 1.22 * 206265.0 * lambda0 * 1E-6 / 0.85 ; beam size (arcsec)
       box = nbeams * beam / pscale                   ; median filter size (pixels)

       print, 'Median filtering'
       mimage = median(image,box)

       newimage = image - mimage   ; should we set negative values to zero?

       sxaddhist, 'Median-filtered by a '+strn(box)+' pixel box.', header
       writefits, datapath+outname[i], newimage, header
       
    endfor

return
end
