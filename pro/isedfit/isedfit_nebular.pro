;+
; NAME:
;   ISEDFIT_NEBULAR()
;
; PURPOSE:
;   Given the number of Lyman-continuum photons and the forbidden
;   emission-line ratios, compute the nebular continuum and
;   emission-line spectrum. 
;
; INPUTS:
;   nlyc - number of Lyman-continuum photons (sec^-1) [NAGE] 
;
; OPTIONAL INPUTS:
;   wave - input/output rest-frame wavelength vector [Angstrom]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   flam_neb - nebular continuum spectrum in units of erg/s/A [NSPEC] 
;
; COMMENTS:
;   See Koleva et al. (2009, Section 3.3) for many details, as well as
;   Leitherer et al. 1999 and Leitherer & Heckman 1995.
;
;   If NLYC is a vector corresponding to different ages then this
;   routine calls itself recursively.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Jul 10, Siena
;
; Copyright (C) 2013, John Moustakas
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

Function isedfit_nebular, nlyc, wave=wave, inst_sigma=inst_sigma, $
  vsigma=vsigma, oiihb=oiihb, oiiihb=oiiihb, niiha=niiha, siiha=siiha

    nage = n_elements(nlyc)
    if nage eq 0L then begin
       doc_library, 'isedfit_nebular'
       return, -1
    endif

; call this routine recursively if there are multiple ages
    if nage gt 1 then begin
       for ii = 0L, nage-1 do begin
          flam_neb1 = isedfit_nebular(nlyc[ii],wave=wave,$
            inst_sigma=inst_sigma,vsigma=vsigma,oiihb=oiihb,$
            oiiihb=oiiihb,niiha=niiha,siiha=siiha)
          if ii eq 0 then flam_neb = flam_neb1 else $
            flam_neb = [[flam_neb],[flam_neb1]]
       endfor
       return, flam_neb
    endif

; some constants    
    light = im_light(/kms)
    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]

; line-width parameters, including instrumental velocity width     
    if n_elements(sigma) eq 0 then sigma = 100.0 ; [km/s]
    if n_elements(inst_sigma) eq 0 then inst_sigma = 70.0 ; [km/s]
    tot_sigma = sqrt(sigma^2+inst_sigma^2) ; [km/s]
    log10sigma = tot_sigma/light/alog(10)  ; line-width [log-10 Angstrom]

; adopt default (typical) forbidden emission-line ratios; I picked
; these values by looking at Fig 1 of Kauffmann+03
    if n_elements(oiiihb) eq 0 then oiiihb = 10.0^(-0.4)
    if n_elements(niiha) eq 0 then niiha = 10.0^(-0.45)
    if n_elements(oiihb) eq 0 then oiihb = 10.0^(0.0)  ; see qa_calibrate_lineratios.pdf
    if n_elements(siiha) eq 0 then siiha = 10.0^(-0.6) ; see qa_calibrate_lineratios.pdf
    
; build the wavelength vector and the output spectrum 
    if n_elements(wave) eq 0 then wave = range(912.0,5E4,5000) ; [Angstrom]
    nwave = n_elements(wave)
    flam_line = wave*0.0
    flam_cont = wave*0.0
    
; gather all the emission-line info we need for the lines of interest,
; starting with the Hydrogen recombination lines; all wavelengths are
; in air!
    lineinfo = {name: '', wave: 0.0, emissivity: 1.0, lum: 0D}

; Hydrogen recombination lines; the emissivities for the first many
; hydrogen transitions are from Storey & Hummer 1995; the file used is
; 'r1b0100.d' where "r" is a prefix, "1" corresponds to "Z=1"
; (hydrogen), "b" corresponds to "case B" (as opposed to case A), and
; "0100" is the nebular temperature (in this case 10,000 K); adopt a
; nebular density of 100 cm^-3; note that the units of the
; emissivities are [erg/s/cm^3]
    nhydrogen = 11
    hydrogen = replicate(lineinfo,nhydrogen)
    hydrogen.name = ['Halpha','Hbeta','Hgamma','Hdelta','Hepsilon',$
      'H6','H7','H8','H9','H10','H11']
    hydrogen.wave = [6562.80,4861.325,4340.464,4101.734,3970.072,$
      3889.049,3835.384,3797.898,3770.630,3750.151,3734.368]
    hydrogen.emissivity = [$
      3.536E-25,$ ; n=3->2 Ha
      1.235E-25,$ ; n=4->2 Hb
      5.784E-26,$ ; n=5->2 Hg
      3.198E-26,$ ; n=6->2 Hd
      1.964E-26,$ ; n=7->2 He
      1.297E-26,$ ; n=8->2 H6
      9.031E-27,$ ; n=9->2 H7
      6.549E-27,$ ; n=10->2 H8
      4.905E-27,$ ; n=11->2 H9
      3.771E-27,$ ; n=12->2 H10
      2.963E-27]  ; n=13->2 H11

; calculate the luminosity of each line (in erg/s/A) given N(Lyc);
; adopted conversions are from Kennicutt 1998 (but see also Leitherer
; & Heckman 1995 and Hao+2011)
    isha = where(hydrogen.name eq 'Halpha',nha)
    ishb = where(hydrogen.name eq 'Hbeta',nhb)
;   niceprint, hydrogen.name, hydrogen[isha].emissivity/hydrogen.emissivity
    for ii = 0, nhydrogen-1 do hydrogen[ii].lum = 1.367D-12*nlyc*$
      hydrogen[ii].emissivity/hydrogen[isha].emissivity/hydrogen[ii].wave
;   struct_print, hydrogen
    
; for the forbidden lines assign empirical line-ratios; adopt average
; doublet line-ratios for [OII] and [SII] corresponding to a density
; of 100/cm^3 and 10^4 K of: <[OII]3726/3729> = 0.75 and
; <[SII]6716/6731>=1.3 (see IM_TEMDEN for details)
    nforbidden = 8
    forbidden = replicate(lineinfo,nforbidden)
;   forbidden.name = ['OII_3727','OIII_5007','NII_6584','SII_6723']
;   forbidden.wave = [3727.4235,5006.843,6583.46,6723.4751]
    forbidden.name = ['OII_3726','OII_3729','OIII_4959','OIII_5007',$
      'NII_6548','NII_6584','SII_6716','SII_6731']
    forbidden.wave = [3726.032,3728.815,4958.911,5006.843,$
      6548.04,6583.46,6716.14,6730.81]

; [OII]    
    oiiratio = 0.75 ; =3726/3729
    is3726 = where(forbidden.name eq 'OII_3726')
    is3729 = where(forbidden.name eq 'OII_3729')
    forbidden[is3726].lum = oiihb*hydrogen[ishb].lum/(1.0+1.0/oiiratio)
    forbidden[is3729].lum = oiihb*hydrogen[ishb].lum/(1.0+oiiratio)

; [SII]    
    siiratio = 1.3 ; =6716/6731
    is6716 = where(forbidden.name eq 'SII_6716')
    is6731 = where(forbidden.name eq 'SII_6731')
    forbidden[is6716].lum = siiha*hydrogen[isha].lum/(1.0+1.0/siiratio)
    forbidden[is6731].lum = siiha*hydrogen[isha].lum/(1.0+siiratio)

; [OIII]
    oiiiratio = 3.0 ; =5007/4959
    is4959 = where(forbidden.name eq 'OIII_4959')
    is5007 = where(forbidden.name eq 'OIII_5007')
    forbidden[is5007].lum = oiiihb*hydrogen[ishb].lum
    forbidden[is4959].lum = oiiihb*hydrogen[ishb].lum/oiiiratio

; [NII]
    niiratio = 3.0 ; =6584/6548
    is6548 = where(forbidden.name eq 'NII_6548')
    is6584 = where(forbidden.name eq 'NII_6584')
    forbidden[is6584].lum = niiha*hydrogen[isha].lum
    forbidden[is6548].lum = niiha*hydrogen[isha].lum/niiratio

; combine the hydrogen and forbidden lines
    line = [hydrogen,forbidden]
    nline = n_elements(line)
;   struct_print, line
    
; build the emission-line spectrum in log-10 wavelength space, since
; each pixel is a fixed width in velocity space (km/s)
    log10pixsz = 0.5D-4
    minlog10wave = alog10(min(line.wave)*0.9)
    nlog10wave = long((alog10(max(line.wave)*1.1)-minlog10wave)/log10pixsz+1)
    log10wave = minlog10wave + dindgen(nlog10wave)*log10pixsz

    lineflux = log10wave*0.0
    for jj = 0, nline-1 do begin
       lineflux += line[jj].lum*exp(-0.5*(log10wave-alog10(line[jj].wave))^2/log10sigma^2)/$
         (sqrt(2.0*!pi)*log10sigma)/alog(10)
    endfor
    flam_line = interpolate(lineflux,findex(10.0^log10wave,wave),missing=0.0)

;   djs_plot, 10^log10wave, lineflux, psym=8, xr=minmax(line.wave)*[0.9,1.1]
;   djs_oplot, wave, flam_line, psym=8, color='green'
    
; scale to a fiducial distance of 10 pc and then add the two
; components together
    flam_line = flam_line/(4.0*!dpi*dist^2)
    flam_cont = flam_cont/(4.0*!dpi*dist^2)
    flam_neb = flam_line + flam_cont
    
;; below are some notes from Pegase and Koleva+09 for computing the
;; nebular *continuum*
;    
;; total recombination coefficient alpha2(T)
;    alpha2 = 2.575D-13 ; Aller+84 [cm^3 s^-1]
;;   alpha2 = 2.616D-13 ; Pegase.2
;
;; assumed densities of He+ and He++ relative to H+, taken from
;; Koleva+09 
;    f_HeI = 0.0897    ; Pegase uses 0.095
;    f_HeII = 1.667D-4 ; Pegase assumes zero
;
;    light = im_light(/ang) ; [Angstrom/s]
;
;; eventually replace this bit of code with my own computations of the
;; total emission coefficient (see for example Koleva+09), but for now
;; just use Pegase's calculations
;;   peg = pegase_hii()
;;   fneb(i)=1.d-40*(g1+g2+HeI*g3+HeII*g4)*c/alphaTe        /lambda(i)**2
;    
;; total emission coefficient
;    gamma_tot = gamma_b + gamma_2p + gamma_HI + $
;      f_HeI*gamma_HeI + f_HeII*gamma_HeII
;
;; see equation (1) in Koleva+09    
;    flam_neb = gamma_tot*nlyc*light/wave^2/alpha2

return, flam_neb
end
