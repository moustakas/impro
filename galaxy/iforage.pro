;+
; NAME:
;       IFORAGE()
;
; PURPOSE:
;       Retrieve useful header information from 1D or 2D spectra. 
;
; INPUTS:
;       speclist - list of FITS spectra
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       forage   - data structure
;
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;       IDL> help, forage('spec2d.fits'), /str
;       IDL> help, forage('spec1d.fits'), /str
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 November 25, U of A - written
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03sep11uofa - replaced FINDFILE() with FILE_SEARCH()
;       jm05jun17uofa - also read the MJD-OBS keyword; improved error
;                       checking; associate POSANGLE not PA as the
;                       slit position angle
;       jm05jul26uofa - bug fix
;       jm07nov29nyu  - got rid of DATAPATH
;       jm09jan19nyu  - added EXPTYPE, ENOISE, and EGAIN header tags
;         (e.g., for Magellan/MIKE data)
;
; Copyright (C) 2003, 2005, 2007, 2009, John Moustakas
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

function readh1d, speclist

    nspec = n_elements(speclist)

    forage = {$
      file:      '', $
      spec2d:    '', $ ; 2D file name
      naxis:     0L, $
      naxis1:    0L, $
      object:    '', $
      type:      '', $
      exptype:   '', $
      galaxy:    '', $
      date:      '', $ ; date observed
      ut:        '', $
      ra:        '', $
      dec:       '', $
      epoch:    0.0, $
      equinox:  0.0, $
      jd:       0.0D,$
      mjd:      0.0D,$
;     pa:       0.0, $ ; slit position angle
      posangle: 0.0, $ ; slit position angle
      airmass:  0.0, $
      zd:       0.0, $ ; zenith distance
      parangle: 0.0, $ ; parallactic angle
      exptime:  0.0, $  ; exposure time
      oexptime: 0.0, $  ; original exposure time
      observat:  '', $
      crval1:   0.0, $
      crpix1:   0.0, $
      cd1_1:    0.0, $
      crval2:   0.0, $
      crpix2:   0.0, $
      cd2_2:    0.0, $
      instrume:  '', $
      telescop:  '', $
      aperture: 0.0, $
      scanlen:  0.0, $
;     extract:  0.0, $
      aperlo:   0.0, $
      aperup:   0.0, $
      aperwid:  0.0, $
      apercen:   0L, $
      medsnr:   0.0, $
      skyshift: 0.0, $
;     xshift:   0.0, $
      zptshift: 0.0, $
      zpterror: 0.0, $
      photflag:  '', $
;     rmaxerr:  0.0, $
;     rmeanerr: 0.0, $
      aminerr:  0.0, $
      amaxerr:  0.0, $
      amederr:  0.0, $
      caliberr: 0.0, $
      abserror: 0.0}

    forage = replicate(forage,nspec)

    tags = tag_names(forage[0])
    ntags = n_elements(tags)
    
    for i = 0L, nspec-1L do begin
       
       h = headfits(speclist[i])
       for j = 0L, ntags-1L do begin
          val = sxpar(h,tags[j],count=count)
          if count ne 0L then forage[i].(j) = val
       endfor

       forage[i].file = speclist[i]
       forage[i].date = strcompress(sxpar(h,'DATE-OBS'),/remove)
       forage[i].type = strcompress(sxpar(h,'IMAGETYP'),/remove)

    endfor

return, forage
end

function readh2d, speclist

    nspec = n_elements(speclist)

    forage = {$
      file:      '', $
      naxis:     0L, $
      naxis1:    0L, $
      naxis2:    0L, $
      object:    '', $
      type:      '', $
      exptype:   '', $
      galaxy:    '', $
      date:      '', $ ; date observed
      ut:        '', $
      ra:        '', $
      dec:       '', $
      epoch:    0.0, $
      jd:       0.0D,$
      mjd:      0.0D,$
;     pa:       0.0, $ ; slit position angle
      posangle: 0.0, $ ; slit position angle
      airmass:  0.0, $
      zd:       0.0, $ ; zenith distance
      parangle: 0.0, $ ; parallactic angle
      exptime:  0.0, $  ; exposure time
      oexptime: 0.0, $  ; original exposure time
      observat:  '', $
      gain:     0.0, $
      rdnoise:  0.0, $
      egain:    0.0, $
      enoise:   0.0, $
      crval1:   0.0, $
      crpix1:   0.0, $
      cd1_1:    0.0, $
      crval2:   0.0, $
      crpix2:   0.0, $
      cd2_2:    0.0, $
      ncombine: 0.0, $
      aperture: 0.0, $
      scanlen:  0.0, $
      zptshift: 0.0, $
      zpterror: 0.0, $
      photflag:  '', $
;     rmaxerr:  0.0, $
;     rmeanerr: 0.0, $
      aminerr:  0.0, $
      amaxerr:  0.0, $
      amederr:  0.0, $
      caliberr: 0.0, $
      abserror: 0.0}

    forage = replicate(forage,nspec)

    tags = tag_names(forage[0])
    ntags = n_elements(tags)
    
    for i = 0L, nspec-1L do begin

       h = headfits(speclist[i])

       for j = 0L, ntags-1L do begin
          val = sxpar(h,tags[j],count=count)
          if (count ne 0L) and (strcompress(val,/remove) ne '') then forage[i].(j) = val
       endfor

       forage[i].file = speclist[i]

       date = sxpar(h,'DATE-OBS',count=ndate)
       if (ndate ne 0L) then forage[i].date = strcompress(date,/remove)
       mjd = sxpar(h,'MJD-OBS',count=nmjd)
       if (nmjd ne 0L) then forage[i].mjd = mjd
       imtype = sxpar(h,'IMAGETYP',count=nimtype)
       if (nimtype ne 0L) then forage[i].type = strcompress(imtype,/remove)

    endfor

return, forage
end

function iforage, speclist

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       doc_library, 'iforage'
       return, -1
    endif

    inspeclist = speclist
;   speclist = file_search(speclist,count=fcount)
;    if (fcount ne n_elements(inspeclist)) then begin ; this needs to be smarter
;       bigindx = lindgen(nspec)
;       match, inspeclist, speclist, suba, subb
;       remove, suba, bigindx
;       splog, 'The following files in SPECLIST were not found:'
;       niceprint, inspeclist[bigindx]
;       return, -1L
;    endif
       
; read the first spectrum: is it 1D or 2D?

    if (file_test(speclist[0],/regular) eq 0L) then begin
       splog, 'File '+speclist[0]+' not found.'
       return, -1
    endif
    
    h = headfits(speclist[0])
    dim = (sxpar(h,'NAXIS'))[0]

    case dim of

       1L: forage = readh1d(speclist)
       2L: forage = readh2d(speclist)

       else: begin
          splog, 'Spectra appear to have more than two dimensions!'
          return, -1
       endelse

    endcase

return, forage
end    
