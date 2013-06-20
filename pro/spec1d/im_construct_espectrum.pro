;+
; NAME:
;   IM_CONSTRUCT_ESPECTRUM()
;
; PURPOSE:
;   Construct a model emission-line spectrum.
;
; CALLING SEQUENCE:
;   espectrum = im_construct_espectrum(wave,sfr=,redshift=,$
;      oiiha=,oiiihb=,niiha=,ebv=,lineres=,linewidth=,$
;      linefit=,/debug,_extra=extra)
;
; INPUTS:
;   wave - *observed* wavelength vector [Angstrom]
;   sfr  - star-formation rate [M_sun/yr]
;
; OPTIONAL INPUTS:
;   redshift  - object redshift (default 0.0)
;   oiiha     - intrinsic [O II] 3727 / H-alpha ratio (1.0)  
;   oiiihb    - intrinsic [O III] 5007 / H-beta ratio (0.8)  
;   niiha     - intrinsic [N II] 6584 / H-alpha ratio (0.5)
;   ebv       - reddening (default 0.0) [mag]
;   lineres   - *observed* FWHM instrumental line-width (default
;               3.0) [Angstrom] 
;   linewidth - intrinsic line-width (default 100) [km/s]
;   extra     - keywords for K_LAMBDA()
;
; KEYWORD PARAMETERS:
;   debug - debug the emission-line spectrum
;
; OUTPUTS:
;   espectrum - final emission-line spectrum [erg/s/cm2/A]; if
;               REDSHIFT=0 then the units are [erg/s/A]
;
; OPTIONAL OUTPUTS:
;   linefit - if DEBUG=1 then return the line-fitting structure
;             results 
;
; PROCEDURES USED:
;   DLUMINOSITY(), RETURN_TBALMER(), MANYGAUSS(), K_LAMBDA(),
;   ILINEFIT(), SPLOG
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Feb 07, U of A - written
;   jm04apr07uofa - extensively tested
;   jm04oct15uofa - added many new inputs; more testing;
;                   documentation added
;   jm05apr13uofa - added EXTRA optional input
;
; Copyright (C) 2004-2005, John Moustakas
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

function gconvolve, x, sigma, _extra=extra
; gaussian convolution
; taken from VDISPFIT

; special case for no smoothing

    if (sigma EQ 0) then return, x

    ksize = round(4*sigma+1) * 2
    xx = findgen(ksize) - ksize/2

    kernel = exp(-xx^2 / (2*sigma^2))
    kernel = kernel / total(kernel)

return, convol(x,kernel,_extra=extra)
end

function onegauss, xval, pp, sigmares=sigmares
; the sigma line-width is comprised of the fixed spectral resolution
; width and the variable intrinsic line width

    sigma_squared = sigmares^2.0 + pp[2]^2.0 ; total Gaussian sigma width [log-Angstrom]
    sigma = sqrt(sigma_squared)

    term1 = exp(-(xval - pp[1])^2 / (2. * sigma_squared ))
    yval = pp[0] * term1 / (sqrt(2.*!pi) * sigma)

;   print, sigmares*3E5*alog(10), pp[2]*3E5*alog(10), sigma*3E5*alog(10)
;   term1 = exp( - (xval - pp[1])^2 / (2. * pp[2]^2) )
;   yval = pp[0] * term1 / (sqrt(2.*!pi) * pp[2])

return, yval
end

function manygauss, xindx, pp, nline=nline, nback=nback, loglam=loglam, $
  sigmares=sigmares, background=background

    yval = 0.0D

    for iline=0, nline-1 do $
      yval = yval + onegauss(loglam[xindx], pp[iline*3:iline*3+2], sigmares=sigmares[iline])
    for iback=0, nback-1 do $
      yval = yval + background[xindx,iback] * pp[nline*3+iback]

;   plot, xindx, yval
;   cc = get_kbrd(1)
    
return, yval
end

function im_construct_espectrum, wave, sfr=sfr, redshift=redshift, oiiha=oiiha, $
  oiiihb=oiiihb, niiha=niiha, ebv=ebv, lineres=lineres, linewidth=linewidth, $
  linefit=linefit, debug=debug, _extra=extra

    light = 2.99792458D5 ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    sfr2energy = 1.0/7.9D-42 ; (erg/s) / (M_sun/yr) [Kennicutt 1998]

    npix = n_elements(wave)
    if (npix eq 0L) then begin
       print, 'Syntax - espectrum = im_construct_espectrum(wave,sfr=,redshift=,$'
       print, '   oiiha=,oiiihb=,niiha=,ebv=,lineres=,linewidth=,linefit=,/debug,$'
       print, '   _extra=extra)'
       return, -1L
    endif
    
    if (n_elements(sfr) eq 0L) then begin
       splog, 'WARNING: SFR is zero.'
       sfr = 0.0
    endif
       
    if (n_elements(redshift) eq 0L) then begin
       redshift = 0.0D 
       fluxarea = 1.0
    endif else begin
       dlum = dluminosity(redshift,/cm)
       fluxarea = 4.0*!dpi*dlum*dlum
    endelse

    if (n_elements(oiiha) eq 0L) then oiiha = 1.0D
    if (n_elements(oiiihb) eq 0L) then oiiihb = 0.8D
    if (n_elements(niiha) eq 0L) then niiha = 0.5D
    if (n_elements(ebv) eq 0L) then ebv = 0.0D
    if (n_elements(lineres) eq 0L) then lineres = 3.0
    if (n_elements(linewidth) eq 0L) then linewidth = 100.0 

; fix the luminosity of H-alpha to the SFR then reduce the Balmer
; lines by the appropriate Balmer decrement

    HaHb = 2.86
    HgHb = 0.468
    HdHb = 0.259
;   HaHb = return_tbalmer(/HaHb)
;   HgHb = return_tbalmer(/HgHb)
;   HdHb = return_tbalmer(/HdHb)
    
    lum_Ha = sfr*sfr2energy
    lum_Hb = lum_Ha/HaHb
    lum_Hg = lum_Hb*HgHb
    lum_Hd = lum_Hb*HdHb
    
; initialize the forbidden emission lines according to the appropriate
; ratios

    oiiratio = 0.8D ; intrinsic 3729/3726 ratio (electron-density dependent)  

    lum_oii_3727 = lum_Ha*oiiha
    lum_oii_3726 = lum_oii_3727/(1+oiiratio)
    lum_oii_3729 = lum_oii_3727*oiiratio/(1+oiiratio)
    
    lum_oiii_5007 = lum_Hb*oiiihb
    lum_oiii_4959 = lum_oiii_5007/3.0

    lum_nii_6584 = lum_Ha*niiha
    lum_nii_6548 = lum_nii_6584/3.0

; initialize the emission lines
    
    linewave = [$
      3726.032D,$ ; [O II]
      3728.815,$ ; [O II]
      4958.911,$ ; [O III]
      5006.843,$ ; [O III]
      4101.734,$ ; Hd
      4340.464,$ ; Hg
      4861.325,$ ; Hb
      6562.800,$ ; Ha
      6548.04, $ ; [N II]
      6583.46]   ; [N II]
    nline = n_elements(linewave)

    area = linewave*0.0 ; this array *must* match LINEWAVE! [erg/s]
    area[0] = lum_oii_3726
    area[1] = lum_oii_3729
    area[2] = lum_oiii_4959
    area[3] = lum_oiii_5007
    area[4] = lum_Hd
    area[5] = lum_Hg
    area[6] = lum_Hb
    area[7] = lum_Ha
    area[8] = lum_nii_6548
    area[9] = lum_nii_6584
    
    amplitude = area/linewave/alog(10.0)   ; kludge for ILINEFIT() [erg/s/A]
    logvwidth = linewidth/light/alog(10.0) ; intrinsic line width [log-km/s]
    
    lineinfo = replicate({value: 0.0D},3*nline)
    lineinfo[0+lindgen(nline)*3].value = amplitude        ; amplitude
    lineinfo[1+lindgen(nline)*3].value = alog10(linewave) ; rest wavelength [log-A]
    lineinfo[2+lindgen(nline)*3].value = logvwidth        ; intrinsic line width [log-km/s]

; ---------------------------------------------------------------------------
; NEW COMMENTS: Construct the emission line spectrum using only the
; intrinsic linewidth (ie, SIGMARES=0); below, using
; IM_GAUSS_BROADEN() we broaden the full spectrum to the *observed*
; (desired) instrumental resolution; also note that we define VSIGRES
; as the *observed* instrumental resolution in [km/s] for use in
; ILINEFIT().  
; ---------------------------------------------------------------------------
; OLD COMMENTS: convert the instrumental resolution (LINERES) to km/s,
; construct the emission-line spectrum; finally redden by E(B-V); NB:
; I am using the *observed* LINERES (rather than RESTLINERES) because
; all we want is a Gaussian of a particular line-width, with no
; reference as to the redshift; instrumental resolution is an
; *observed* quantity
; ---------------------------------------------------------------------------

    redlinewave = linewave*(1+redshift)
    vsigres = light*lineres/redlinewave/fwhm2sig ; *observed* instrumental resolution [km/s]
;   sigmares = vsigres/light/alog(10.0)          ; instrumental resolution [log-km/s]
    sigmares = linewave*0.0

    restwave = wave/(1+redshift)
    logrestwave = alog10(restwave)
    xindx = lindgen(npix)

    espectrum1 = manygauss(xindx,lineinfo.value,nline=nline,$ ; [erg/s/A]
      nback=0,loglam=logrestwave,sigmares=sigmares)

; redden, redshift, then finally broaden to the desired *observed*
; spectral resolution
    
    espectrum2 = espectrum1*10^(-0.4*ebv*k_lambda(restwave,/odonnell))
;   espectrum3 = espectrum2/fluxarea              ; [erg/s/cm2/A]
    espectrum3 = espectrum2/fluxarea/(1+redshift) ; [erg/s/cm2/A]
    espectrum = im_gauss_broaden(wave,espectrum3,0.0,lineres)

    if keyword_set(debug) then begin

       restespec = espectrum*(1+redshift)
       norm = max(restespec)

       
       
       linefit = ilinefit(restwave,restespec/norm,espectrum*0.0+0.1,$
         linewave,vsigres,$
         ,specfit=specfit,zguess=0.0)
;      linefit = ilinefit(espectrum/norm,wave,linewave,vsigres,$
;    invvar=espectrum*0.0+0.1,specfit=specfit,zguess=redshift)

       splog, '[O II] ratio =  ', (linefit[0].linearea+linefit[1].linearea)/linefit[7].linearea, oiiha
       splog, '[O III] ratio = ', linefit[3].linearea/linefit[2].linearea, 3.0
       splog, 'Ha/Hb =         ', linefit[7].linearea/linefit[6].linearea, hahb
       splog, 'Hb/Hg =         ', linefit[6].linearea/linefit[5].linearea, 1.0/hghb, $
         get_ebv(linefit[6].linearea/linefit[5].linearea,/hbeta)
       splog, 'Hb/Hd =         ', linefit[6].linearea/linefit[4].linearea, 1.0/hdhb

       window, 0, xs=800, ys=450
       plot, wave, espectrum, ps=10, xsty=3, ysty=3, thick=1.0, $
         xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, $
         xtitle='Wavelength', ytitle='Flux'
       djs_oplot, wave, norm*specfit, ps=10, color='yellow', thick=2.0

       linefit.linearea = linefit.linearea*norm
       linefit.linearea_err = linefit.linearea_err*norm
       linefit.linebox = linefit.linebox*norm
       linefit.linebox_err = linefit.linebox_err*norm
       linefit.linecontlevel = linefit.linecontlevel*norm
       linefit.linecontlevel_err = linefit.linecontlevel_err*norm

    endif
    
return, espectrum
end
