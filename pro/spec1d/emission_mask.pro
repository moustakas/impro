;+
; NAME:
;       EMISSION_MASK()
;
; PURPOSE:
;       Mask optical nebular emission lines.
;
; INPUTS:
;       wave  - wavelength vector [Angstrom]
;
; OPTIONAL INPUTS:
;       z         - dimensionless redshift
;       width     - width of the mask around each emission line
;                   [Angstrom] (default 15.0)
;       bluewidth - width of the mask around each blue emission line  
;                   [Angstrom] (default 5.0); only used if BLUEMASK=1 
;       abswidth  - width of the mask around strong absorption lines
;                   [Angstrom] (default 10.0); only used if ABSMASK=1 
;       spectrum  - spectrum corresponding to WAVE
;       ivarspectrum - inverse variance spectrum corresponding to WAVE
;       linelist  - *rest* wavelengths of nebular lines to mask 
;       sigrej_cosmic - rejection threshold if COSMIC=1
;       extra     - keywords for DJS_ITERSTAT
;
; KEYWORD PARAMETERS:
;       telluric  - flag wavelengths corresponding to telluric
;                   absorption lines in MASK
;       cosmic    - attempt to identify and flag cosmic rays (deviant
;                   pixels) in SPECTRUM and add them to MASK 
;       negative  - mask pixels that are negative
;       bluemask  - mask the high-order Balmer lines in the blue part
;                   of the optical spectrum
;       absmask   - mask strong absorption lines
;       noskymask - do not mask sky wavelengths 
;       qsomask   - mask QSO lines
;       gauss     - fit a Gaussian to each nebular line and mask
;                   within 3-sigma (overides WIDTH and requires
;                   SPECTRUM and IVARSPECTRUM)
;
; OUTPUTS:
;       mask  - pixel mask for WAVE (1B is good and 0B is bad)
;
; OPTIONAL OUTPUTS:
;       good  - indices in WAVE that have not been masked
;       bad   - indices in WAVE that have been masked
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 April 8, U of A
;       jm03nov27uofa - made TELLURIC_MASK() its own function 
;       jm04jan05uofa - added NEGATIVE keyword
;       jm04jan14uofa - also mask sky wavelengths
;       jm04may13uofa - added BLUEMASK keyword
;       jm05jan28uofa - added BLUEWIDTH optional input; mask the blue
;                       lines by a smaller amount relative to the
;                       other lines
;       jm05aug02uofa - added NOSKYMASK keyword
;       jm06feb08uofa - added SIGREJ_COSMIC optional input
;       jm06feb16uofa - added QSOMASK keyword
;       jm07mar02nyu  - added QMASK optional output and added more QSO
;                       emission-line wavelengths 
;       jm07sep29nyu  - added ABSMASK and ABSWIDTH keywords 
;       jm09jan24nyu  - added LINELIST and IVARSPECTRUM optional
;                       inputs and GAUSS keyword 
;
; Copyright (C) 2002-2005, 2007, 2009, John Moustakas
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

function emission_mask, wave, z=z, width=width, bluewidth=bluewidth, abswidth=abswidth, $
  qsowidth=qsowidth, spectrum=spectrum, ivarspectrum=ivarspectrum, linelist=linelist1, $
  good=good, bad=bad, sigrej_cosmic=sigrej_cosmic, telluric=telluric, cosmic=cosmic, $
  negative=negative, gauss=gauss, bluemask=bluemask, absmask=absmask, noskymask=noskymask, $
  qsomask=qsomask, qmask=qmask, crmask=crmask, skymask=skymask, tellmask=tellmask, _extra=extra

    nwave = n_elements(wave)
    if nwave eq 0L then begin
       doc_library, 'emission_mask'
       return, -1
    endif

    nspec = n_elements(spectrum)

    if n_elements(z) eq 0L then z = 0.0
    if n_elements(width) eq 0L then width = 15.0        ; mask width [Angstrom]
    if n_elements(bluewidth) eq 0L then bluewidth = 5.0 ; blue mask width [Angstrom]
    if n_elements(abswidth) eq 0L then abswidth = 20.0  ; absorption-line mask width [Angstrom]
    if n_elements(qsowidth) eq 0L then qsowidth = 100.0 ; QSO mask width [Angstrom]
    
    mask = make_array(nwave,/byte,value=1)

    if keyword_set(gauss) and ((n_elements(spectrum) eq 0L) or $
      (n_elements(ivarspectrum) eq 0L)) then begin
       splog, 'GAUSS keyword requires SPECTRUM and IVARSPECTRUM input'
       return, mask
    endif

    if (n_elements(linelist1) eq 0L) then begin
       linelist1 = [$
         3425.868,$                ; OII
         3726.032,$                ; OII
         3728.815,$                ; OII
         4101.734,$                ; H-delta
         4340.46 ,$                ; H-gamma
         4363.21 ,$                ; OIII
         4861.33 ,$                ; H-beta
         4958.91 ,$                ; OIII
         5006.84 ,$                ; OIII
         5875.96 ,$                ; HeI
         5890.0  ,$                ; Na D doublet
         5896.0  ,$                ; Na D doublet
         6300.30 ,$                ; OI
         6548.04 ,$                ; NII
         6562.80 ,$                ; H-alpha
         6583.46 ,$                ; NII
         6716.14,$                 ; SII
         6730.81$                  ; SII
       ]
    endif 
    linelist = (1+z)*linelist1 
       
    for i = 0L, n_elements(linelist)-1L do begin
; note: the algorithm will not properly mask emission lines on the
; edges because I require that the line be fully contained in the
; spectrum 
       mflag = 0
       if keyword_set(gauss) then begin
          zoom_width = 5.0
          linepix = where((wave gt (linelist[i]-zoom_width)) and $
            (wave lt (linelist[i]+zoom_width)),npix)
          if (npix ne 0L) then begin
             estimates = [im_max(spectrum[linepix],sigrej=3.0),$
               linelist[i],1.0,im_median(spectrum[linepix],sigrej=3.0)]
             yfit = mpfitpeak(wave[linepix],spectrum[linepix],params,$
               /positive,/gaussian,perror=perror,estimates=estimates,$
               weights=ivarspectrum[linepix])
             if (params[0]/(perror[0]+(perror[0] eq 0.0))*(perror[0] ne 0.0) gt 3.0) then begin
;               djs_plot, wave[linepix],spectrum[linepix]
;               djs_oplot, wave[linepix], yfit, color='red'            
                mask = mask and ((wave lt params[1]-3.0*params[2]) or $
                  (wave gt params[1]+3.0*params[2]))
             endif else mflag = 1
          endif else mflag = 1
       endif else mflag = 1
       if mflag then mask = mask and ((wave lt linelist[i]-width) or $
         (wave gt linelist[i]+width))
    endfor 

    skymask = make_array(nwave,/byte,value=1)
    skylist = [5577.339,5889.950,6300.32,6363.81] ; sky wavelengths

    for i = 0L, n_elements(skylist)-1L do $
      skymask = skymask and ((wave lt skylist[i]-width) or (wave gt skylist[i]+width))

    if (not keyword_set(noskymask)) then mask = (mask + skymask) gt 1B
    
; mask blue pixels    
    
    if keyword_set(bluemask) then begin

       bmask = make_array(nwave,/byte,value=1)

       bluelinelist = (1+z)*[$
         3734.368,$             ; H11
         3750.151,$             ; H10
         3770.630,$             ; H9
         3797.898,$             ; H8
         3835.384,$             ; H7
         3869.06, $             ; NeIII
         3889.049,$             ; H6
         3970.072 $             ; H-epsilon
         ]

       for i = 0L, n_elements(bluelinelist)-1L do $
         bmask = bmask and ((wave lt bluelinelist[i]-bluewidth) or $
         (wave gt bluelinelist[i]+bluewidth))
       mask = (mask + bmask) eq 2B
       
    endif
    
; mask strong absorption lines
    
    if keyword_set(absmask) then begin

       amask = make_array(nwave,/byte,value=1)

       abslinelist = (1+z)*[$
         3930.0,$ ; Ca
         3968.0$  ; Ca
         ]

       for i = 0L, n_elements(abslinelist)-1L do $
         amask = amask and ((wave lt abslinelist[i]-abswidth) or $
         (wave gt abslinelist[i]+abswidth))
       mask = (mask + amask) eq 2B
       
    endif
    
; mask QSO lines
    
    if keyword_set(qsomask) then begin

       qmask = make_array(nwave,/byte,value=1)

       qsolinelist = (1+z)*[$
         1218.0,$
         1305.0,$
         1400.0,$
         1549.0,$
         1909.0,$
         2799.495$
         ]

       for i = 0L, n_elements(qsolinelist)-1L do $
         qmask = qmask and ((wave lt qsolinelist[i]-qsowidth) or $
         (wave gt qsolinelist[i]+qsowidth))
       mask = (mask + qmask) eq 2B
       
    endif
    
; flag telluric features
    
    if keyword_set(telluric) then begin
       tellmask = telluric_mask(wave)
       mask = (mask + tellmask) eq 2B
    endif else begin
       if arg_present(tellmask) then tellmask = mask*0B+1B
    endelse

; attempt to remove cosmic-rays if requested    

    if keyword_set(cosmic) and (nspec ne 0L) then begin

       if (nspec ne nwave) then begin
          splog, 'SPECTRUM and WAVE must have the same number of elements.'
       endif else begin

;         res = djs_avsigclip(reform(spectrum,nspec,1),0,inmask=(mask eq 0B),$
;           outmask=crmask,sigrej=3.0)
;         crmask = crmask eq 0B

          if (n_elements(sigrej_cosmic) eq 0L) then sigrej_cosmic = 3.0

          noelines = where(mask eq 1B)
          djs_iterstat, spectrum[noelines], sigrej=sigrej_cosmic, mask=crmask1, _extra=extra

;         crmask1 = (smooth(float(crmask1 eq 0),2) gt 0.0) eq 0B ; grow the mask

          mask[noelines] = (mask[noelines] + crmask1) gt 1B

          crmask = make_array(nwave,/byte,value=1)
          crmask[noelines] = crmask1
          
       endelse
    endif

    if keyword_set(negative) and (nspec ne 0L) then begin

       if (nspec ne nwave) then begin
          splog, 'SPECTRUM and WAVE must have the same number of elements.'
       endif else begin

          neg = where(spectrum lt 0.0,nneg)
          if (nneg ne 0L) then mask[neg] = 0B

       endelse
       
    endif

    good = where(mask eq 1B,ngood,comp=bad,ncomp=nbad)

return, mask
end
