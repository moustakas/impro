;+
; NAME:
;       IMATCH_TEMPLATES()
;
; PURPOSE:
;       Match a set of templates to a data spectrum in terms of
;       wavelength spacing and spectral resolution.
;
; CALLING SEQUENCE:
;          newflux = imatch_templates(starflux,starwave,$
;             newwave,starres=,newres=)
;
; INPUTS:
;       starflux - stellar template flux array [NSTARFLUX,NSTAR] 
;       starwave - wavelength vector corresponding to STARFLUX
;                  [NSTARFLUX]
;       newwave  - desired output wavelength spacing (i.e., of the
;                  data spectrum)
;
; OPTIONAL INPUTS:
;       starres  - FWHM spectral resolution of STARFLUX [Angstrom] 
;       newres   - FWHM spectral resolution of NEWFLUX [Angstrom]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       newflux  - resampled and broadened STARFLUX
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;       The templates will be padded with zeros where there is no
;       wavelength overlap.  
;
;       This routine currently only supports a single Gaussian sigma
;       as the broadening kernel.  The slow-down is significant for a
;       wavelength-varying kernel. 
;
; PROCEDURES USED:
;       COMBINE1FIBER, IM_GAUSS_BROADEN()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 8, U of A
;       jm03jun6uofa - excised from IFITSPEC()
;       jm03sep15uofa - bug fix on the broadening resolution thanks to
;                       C. Tremonti
;       jm04jan14uofa - if STARRES is negative then do not convolve
;       jm04feb04uofa - broaden *then* interpolate
;       jm04mar02uofa - allow COMBINE1FIBER to determine the inverse
;                       variance from the data themselves
;       jm04may09uofa - crop the wavelength range before broadening to
;                       gain speed
;
; Copyright (C) 2003-2004, John Moustakas
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

function imatch_templates, starflux, starwave, newwave, $
  starres=starres, newres=newres

    npix = n_elements(newwave)

    if (n_elements(starflux) eq 0L) or (n_elements(starwave) eq 0L) or (npix eq 0L) then begin
       print, 'Syntax - newflux = imatch_templates(starflux,starwave,$'
       print, '   newwave,starres=,newres=)'
       return, -1L
    endif
    
    sdim = size(starflux,/n_dimension)
    if sdim ne 2L then begin
       print, 'STARFLUX must have two dimensions.'
       return, -1L
    endif
    
    ssize = size(starflux,/dimension)
    nstar = ssize[1]
    nstarflux = ssize[0]

    if n_elements(starwave) ne nstarflux then begin
       print, 'STARFLUX and STARWAVE must have the same number of pixels.'
       return, -1L
    endif

    nstarres = n_elements(starres)
    
; broaden unless STARRES is negative

    if (starres[0] gt 0.0) then begin

; to expedite the broadening crop STARWAVE and STARFLUX to be within
; 10 pixels of the edges of NEWWAVE

       linterp, newwave, newres, starwave, newstarres, /nointerp

       broadflux = starflux
       get_element, starwave, minmax(newwave), ww
       w1 = (ww[0]-10L)>0L
       w2 = (ww[1]+10L)<(nstarflux-1L)

       if (nstarres ne 1L) then instarres = starres[w1:w2] else instarres = starres
       
       broadflux[w1:w2,*] = im_gauss_broaden(starwave[w1:w2],starflux[w1:w2,*], $
         instarres,newstarres[w1:w2])

    endif else broadflux = starflux

; next re-sample the stellar templates onto the new wavelength spacing 

    newflux = fltarr(npix,nstar)

    for istar = 0L, nstar-1L do begin

       combine1fiber, alog10(starwave), double(broadflux[*,istar]), $
         newloglam=alog10(newwave), newflux=outflux
       newflux[*,istar] = float(outflux)

    endfor

return, newflux
end

