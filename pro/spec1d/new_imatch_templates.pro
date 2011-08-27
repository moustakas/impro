;+
; NAME:
;       NEW_IMATCH_TEMPLATES()
;
; PURPOSE:
;       Match a set of templates to a galaxy spectrum in terms of
;       wavelength spacing, instrumental resolution, and velocity
;       dispersion.
;
; CALLING SEQUENCE:
;          new_templateflux = new_imatch_templates(templateflux,$
;             templatewave,galaxywave,templateres=,galaxyres=,$
;             galaxyvdisp=)
;
; INPUTS:
;       templateflux - template flux [NTEMPLATEPIX,NTEMPLATE] 
;       templatewave - corresponding rest wavelength [Angstrom, 
;                      NTEMPLATEPIX]; assume *constant* pixel size in
;                      A/pixel 
;       galaxywave   - galaxy rest wavelength [Angstrom, NPIX]
;
; OPTIONAL INPUTS:
;       templateres  - *scalar* FWHM spectral resolution of
;                      TEMPLATEFLUX [Angstrom]
;       galaxyres    - desired output, *scalar*, rest FWHM
;                      instrumental resolution of NEWFLUX [Angstrom],
;                      which must be larger than TEMPLATERES
;       galaxyvdisp  - *scalar*, rest-frame galaxy velocity dispersion 
;                      [km/s]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       new_templateflux - resampled and broadened TEMPLATEFLUX
;                          [NPIX,NTEMPLATE] 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The templates will be padded with zeros where there is no
;       wavelength overlap.  
;
;       This routine is optimized/written for the following
;       circumstances: scalar template resolution, scalar galaxy
;       instrumental resolution, scalar galaxy velocity
;       dispersion, and constant template pixel size.  First we
;       degrade the templates to the instrumental resolution of the
;       data in wavelength space, and then resample to constant
;       log-lambda and ... 
;
; PROCEDURES USED:
;       COMBINE1FIBER
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 May 09, U of A - re-written version of
;          IMATCH_TEMPLATES() 
;
; Copyright (C) 2004, John Moustakas
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

function new_imatch_templates, templateflux, templatewave, galaxywave, $
  templateres=templateres, galaxyres=galaxyres, galaxyvdisp=galaxyvdisp

    light = 2.99792458D5               ; speed of light [km/s]
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0)) ; = 2.35

    npix = n_elements(galaxywave)

    if (n_elements(templateflux) eq 0L) or (n_elements(templatewave) eq 0L) or (npix eq 0L) then begin
       print, 'Syntax - new_templateflux = new_imatch_templates(templateflux'
       print, '   templatewave,galaxywave,templateres=,galaxyres=,galaxyvdisp=)'
       return, -1L
    endif
    
    sdim = size(templateflux,/n_dimension)
    if sdim ne 2L then begin
       splog, 'TEMPLATEFLUX must have two dimensions.'
       return, -1L
    endif
    
    ssize = size(templateflux,/dimension)
    ntemplate = ssize[1]
    ntemplatepix = ssize[0]

    if n_elements(templatewave) ne ntemplatepix then begin
       splog, 'TEMPLATEFLUX and TEMPLATEWAVE must have the same number of pixels.'
       return, -1L
    endif

    if (n_elements(templateres) ne 1L) and (n_elements(galaxyres) ne 1L) $
      and (n_elements(galaxyvdisp) ne 1L) then begin
       splog, 'TEMPLATERES, GALAXYRES, and GALAXYVDISP all must be scalars.'
       return, -1L
    endif

    if (templateres gt galaxyres) then begin
       splog, 'TEMPLATERES must be smaller than GALAXYRES.'
       return, -1L
    endif

; to speed up the broadening crop TEMPLATEWAVE and TEMPLATEFLUX to be
; within 10 pixels of the wavelength range of GALAXYWAVE

    cropflux = templateflux
    get_element, templatewave, minmax(galaxywave), ww
    w1 = (ww[0]-10L)>0L
    w2 = (ww[1]+10L)<(ntemplatepix-1L)

    crop_templatewave = templatewave[w1:w2]
    crop_templateflux = templateflux[w1:w2,*]
    
; assume the template pixel size is constant; choose a constant
; velocity pixel size of ~30 km/s (1D-4 in log-lambda)

    template_pixsize = crop_templatewave[1]-crop_templatewave[0]
    sigma_instr = sqrt(galaxyres^2.0 - templateres^2.0)/fwhm2sig/template_pixsize ; [pixels]

    coeff0 = alog10(min(crop_templatewave))
    coeff1 = 1D-4

    lognpix = long(1.0D + (alog10(max(crop_templatewave))-coeff0)/coeff1)
    newlogwave = coeff0 + coeff1*findgen(lognpix)
    newwave = 10.0^newlogwave

    sigma_vdisp = galaxyvdisp/coeff1/light ; [pixels]

; loop on each template

    new_templateflux = galaxywave # replicate(1.0,ntemplate)

    for j = 0L, ntemplate-1L do begin

; first, broaden FROM TEMPLATERES to GALAXYRES

       crop_templateflux[*,j] = gconvolve(reform(crop_templateflux[*,j]),$
         sigma_instr,/edge_truncate)

; next, resample the templates to constant pixel size in km/s and
; broaden to GALAXYVDISP

       combine1fiber, alog10(crop_templatewave), reform(crop_templateflux[*,j]), $
         newloglam=newlogwave, newflux=vbroadflux
       vbroadflux = gconvolve(vbroadflux,sigma_vdisp,/edge_truncate)

; finally, resample onto the pixel size and wavelength spacing of the
; galaxy (GALAXYWAVE)

       combine1fiber, newlogwave, vbroadflux, newloglam=alog10(galaxywave), $
         newflux=newflux
       new_templateflux[*,j] = newflux
       
    endfor

    plot, templatewave, templateflux[*,3], ps=10, xsty=3, ysty=3, $
      xrange=[4000,5100]
    djs_oplot, galaxywave, new_templateflux[*,3], ps=10, color='red', thick=3.0

return, new_templateflux
end

