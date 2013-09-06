;+
; NAME:
;	IM_FITS_CUBE()
;
; PURPOSE:
;	Rad many FITS images into a structure cube.
;
; CALLING SEQUENCE:
;	cube = fits_cube(flist, [path=])
;
; INPUTS:
;	flist - string array of FITS file names
;
; OPTIONAL INPUTS:
;	path  - the data path to the FITS files
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	cube  - structure cube with an image and header field 
;
; OPTIONAL OUTPUTS:
;	head  - header of the first image
;
; COMMENTS:
;	Extra keywords can be passed to READFITS.  This routine is
;	based on DJS_READMANY, but has been expanded to return the
;	FITS header as a structure field.
;
; PROCEDURES USED:
;	READFITS()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Aug 12, U of A, written
;       jm02nov1uofa - added additional error checking
;       jm02nov27uofa - generalized to read in one-dimensional FITS
;                       files
;       jm04nov05uofa - updated the HEADER string array from 500 to
;                       1000 elements
;       jm13may27siena - use GZ_MRDFITS in lieu of READFITS()
;
; Copyright (C) 2001-2002, 2004, 2013, John Moustakas
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

function im_fits_cube, flist, path=path, _extra=extra

    nflist = n_elements(flist)
    if nflist eq 0 then begin
       doc_library, 'im_fits_cube'
       return, -1
    endif

    if (size(flist[0],/type) ne 7) or (strmatch(flist[0],'*.fits*') eq 0) then begin
       print, 'File list must be type string FITS files.'
       return, -1
    endif

    if n_elements(path) eq 0 then path = get_pwd()

    pushd, path
    
; read in the first image

    if file_test(flist[0]+'*',/regular) eq 0L then begin
       splog, 'FITS file '+flist[0]+' not found.'
       return, -1L
    endif

    image = gz_mrdfits(flist[0],0,head,_extra=extra,/silent)
    ndim = size(image,/n_dim)
    imsz = size(image,/dim)
    
    if (ndim ne 1L) and (ndim ne 2L) then begin
       splog, 'FITS file must be one- or two-dimensional!'
       return, -1L
    endif
    
; create the data cube, assuming all the images are the same size.
; store the headers in an array of pointers

    if ndim eq 1L then $
      template = {fname: '', image: make_array(imsz[0],/float), nhead: 0, $
      header: strarr(1000)} else $
      template = {fname: '', image: make_array(imsz[0],imsz[1],/float), nhead: 0, $
      header: strarr(1000)}
    cube = replicate(template,nflist)

    hgood = where(strcompress(head,/remove) ne '',nhgood)
    cube[0].fname = flist[0]
    cube[0].image = image
    cube[0].nhead = nhgood
    cube[0].header[0:nhgood-1] = head[hgood]

    for i = 1, nflist-1L do begin

       if file_test(flist[i],/regular) eq 0L then begin
          splog, 'FITS file '+flist[i]+' not found.'
          return, -1L
       endif
    
       image = gz_mrdfits(flist[i],0,head,_extra=extra,/silent)

       szt = size(image,/dim)
       sdim = size(image,/n_dim)

       if (sdim ne ndim) then begin
          splog, 'FITS files have incompatible dimensions!'
          return, -1L
       endif

       if sdim eq 1L then begin
          if (imsz[0] ne szt[0]) then begin
             splog, 'FITS files are not the same size!'
             return, -1L
          endif
       endif else begin
          if ((imsz[0] ne szt[0]) or (imsz[1] ne szt[1])) then begin
             splog, 'FITS files are not the same size!'
             return, -1L
          endif
       endelse 

       hgood = where(strcompress(head,/remove) ne '',nhgood)
       cube[i].fname = flist[i]
       cube[i].image = image
       cube[i].nhead = nhgood
       cube[i].header[0:nhgood-1] = head[hgood]
    endfor

    popd
    
return, cube
end
