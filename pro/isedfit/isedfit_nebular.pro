;+
; NAME:
;   ISEDFIT_NEBULAR()
;
; PURPOSE:
;   Compute the nebular continuum and emission-line spectrum given the
;   number of Lyman-continuum photons and the forbidden line-ratios
;   calibrated using ISEDFIT_CALIBRATE_LINERATIOS.
;
; INPUTS:
;   nlyc - number of Lyman-continuum photons (sec^-1) [NSPEC]
;
; OPTIONAL INPUTS:
;   wave - input/output rest-frame wavelength vector [NPIX, Angstrom];
;     if WAVE is given then you need to ensure adequate sampling
;     around the wavelengths of the lines!
;   inst_vsigma - instrumental velocity dispersion (default 0 km/s)
;   vsigma - intrinsic (kinematic) velocity dispersion (default 100 km/s) 
;   oiiihb - the forbidden emission-line strengths are calibrated
;     empirically (see ISEDFIT_CALIBRATE_LINERATIOS) to the
;     [OIII]/H-beta line-ratio; the default value is 10.0^(-0.4
;     dex)=0.398, which is typical of a galaxy with a ~solar
;     (gas-phase) metallicity [scalar or NSPEC vector]
;
; KEYWORD PARAMETERS:
;   nospectrum - do not return the spectrum (used by
;     ISEDFIT_BUILD_MONTEGRIDS to generate an emission-line friendly
;     wavelength array)
;
; OUTPUTS:
;   flam_neb - nebular continuum spectrum in units of erg/s/cm^2/A
;     placed at a fiducial distance of 10 pc [NPIX,NSPEC]
;
; OPTIONAL OUTPUTS:
;   line - data structure containing all the details of the
;     emission-line spectrum [NLINE]
;   flam_line - pure nebular emission-line spectrum
;     [NPIX,NSPEC,erg/s/cm^2/A] 
;   flam_cont - pure nebular continuum spectrum
;     [NPIX,NSPEC,erg/s/cm^2/A]  
;
; COMMENTS:
;   If NLYC is a vector corresponding to different desired output
;   spectra then this routine calls itself recursively.  In that case,
;   OIIIHB can be a scalar or a vector.
;
;   Note that FLAM_NEB = FLAM_LINE + FLAM_CONT
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Aug 01, Siena
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

Function isedfit_nebular, nlyc, wave=wave, inst_vsigma=inst_vsigma, $
  vsigma=vsigma, oiiihb=oiiihb1, line=line, flam_line=flam_line, $
  flam_cont=flam_cont, nospectrum=nospectrum

    common isedfit_neb, branch, nebcont, hdata, fdata
    
    nspec = n_elements(nlyc)
    if nspec eq 0L then begin
       doc_library, 'isedfit_nebular'
       return, -1
    endif

; some constants    
    nHeI = 0.0897               ; assumed n(HeI)/n(HI) abundance ratio (see Kotulla+09)
    dist = 10.0*3.085678D18     ; fiducial distance [10 pc in cm]
    light = im_light(/kms)      ; [km/s]

; line-width parameters, including instrumental velocity width     
    if n_elements(vsigma) eq 0 then vsigma = 100.0          ; [km/s]
    if n_elements(inst_vsigma) eq 0 then inst_vsigma = 70.0 ; approximately 3 A FWHM @5500 A [km/s]
    tot_vsigma = sqrt(vsigma^2+inst_vsigma^2)               ; [km/s]
    log10sigma = tot_vsigma/light/alog(10)                  ; line-width [log-10 Angstrom]

; adopt a fiducial [OIII]/Hbeta line-ratio
    noiiihb = n_elements(oiiihb1)
    if noiiihb ne 0 then begin
       if noiiihb eq 1 then begin
          if nspec eq 1 then oiiihb = oiiihb1 else $
            oiiihb = replicate(oiiihb1,nspec)
       endif else begin
          if noiiihb ne nspec then begin
             splog, 'Dimensions of OIIIHB and NLYC must match!'
             return, -1
          endif
       endelse
    endif else begin
       if nspec eq 1 then oiiihb = 10.0^(-0.4) else $
         oiiihb = replicate(10.0^(-0.4),nspec)
    endelse 

; call this routine recursively
    if nspec gt 1 then begin
       for ii = 0L, nspec-1 do begin
          flam_neb1 = isedfit_nebular(nlyc[ii],wave=wave,line=line,$
            inst_vsigma=inst_vsigma,vsigma=vsigma,oiiihb=oiiihb[ii],$
            nospectrum=nospectrum,flam_line=flam_line1,flam_cont=flam_cont1)
          if ii eq 0 then begin
             flam_neb = flam_neb1 
             flam_line = flam_line1 
             flam_cont = flam_cont1 
          endif else begin
             flam_neb = [[flam_neb],[flam_neb1]]
             flam_line = [[flam_line],[flam_line1]]
             flam_cont = [[flam_cont],[flam_cont1]]
          endelse
       endfor
       return, flam_neb
    endif

; read the branching ratios we'll need to build the doublets
; below; assume default density and temperature (100/cm^3 and 10^4 K) 
    if n_elements(branch) eq 0 then branch = im_branch_ratios()
    
; read the nebular continuum file (see ISEDFIT_CALIBRATE_LINERATIOS)  
    if n_elements(nebcont) eq 0 then begin
       nebfile = getenv('IMPRO_DIR')+'/etc/isedfit_nebular_continuum.fits.gz'
       if file_test(nebfile) eq 0 then begin
          splog, 'Nebular continuum file '+nebfile+' not found!'
          return, -1
       endif
       nebcont = mrdfits(nebfile,1)
    endif

; read the hydrogen and helium emissivities
    if n_elements(hdata) eq 0 then begin
       hfile = getenv('IMPRO_DIR')+'/etc/hydrogen_helium_emissivities.dat'
       if file_test(hfile) eq 0 then begin
          splog, 'Hydrogen/helium emissivities file '+hfile+' not found!'
          return, -1
       endif
       hdata = rsex(hfile)
    endif
    
; read the forbidden emission-line ratios (see
; ISEDFIT_CALIBRATE_LINERATIOS) 
    if n_elements(fdata) eq 0 then begin
       ffile = getenv('IMPRO_DIR')+'/etc/isedfit_forbidden_lineratios.fits.gz'
       if file_test(ffile) eq 0 then begin
          splog, 'Forbidden line-ratios file '+ffile+' not found!'
          return, -1
       endif
       fdata = mrdfits(ffile,1)
    endif

    nhline = n_elements(hdata)
    nfline = n_elements(fdata)
    nline = nhline + nfline

; build the emission-line luminosities, starting with the hydrogen and
; helium lines; calculate the luminosity of each line (in erg/s/A)
; given N(Lyc); adopted conversions are from Kennicutt 1998 (but see
; also Leitherer & Heckman 1995 and Hao+2011)
    hlineinfo = replicate({id: 0, name: '', wave: 0D, flux: 0D, ratio: 0.0},nhline)
    hlineinfo.name = hdata.name
    hlineinfo.wave = hdata.wave
    
    isha = where(strtrim(hdata.name,2) eq 'Halpha',nha)
    ishb = where(strtrim(hdata.name,2) eq 'Hbeta',nhb)
    if nha eq 0 or nhb eq 0 then message, 'Problem with the hydrogen/helium data file!'
;   niceprint, hdata.name, hdata[isha].emissivity/hdata.emissivity

    hlineinfo.ratio = hdata.emissivity/hdata[ishb].emissivity
    for ii = 0, nhline-1 do begin
       if strmatch(strtrim(hlineinfo[ii].name,2),'HeI_*') then $
         factor = nHeI else factor = 1.0
       hlineinfo[ii].flux = factor*1.367D-12*nlyc[0]*$ ; [erg/s/A]
         hdata[ii].emissivity/hdata[isha].emissivity/hdata[ii].wave
    endfor
;   struct_print, hlineinfo

; now the forbidden lines; see ISEDFIT_CALIBRATE_LINERATIOS for
; details, but basically *all* the forbidden lines are tied to the
; [OIII]/Hbeta ratio
    flineinfo = replicate({id: 0, name: '', wave: 0D, flux: 0D, ratio: 0.0},nfline)
    flineinfo.name = strtrim(fdata.name,2)
    flineinfo.wave = fdata.wave
    
    for ii = 0, nfline-1 do begin
       ratio = poly(alog10(oiiihb[0]),fdata[ii].coeff)+randomn(seed)*fdata[ii].scatter
       flineinfo[ii].ratio = 10D^ratio
       flineinfo[ii].flux = hlineinfo[ishb].flux*flineinfo[ii].ratio
    endfor

; add [OIII] 5007, with no scatter
    oiii = im_empty_structure(flineinfo[0],ncopies=1)
    oiii.name = '[OIII]_5007'
    oiii.wave = 5006.842D
    oiii.ratio = oiiihb[0]
    oiii.flux = hlineinfo[ishb].flux*oiiihb[0]
    flineinfo = [flineinfo,oiii]

; now add the doublets; adopt average doublet line-ratios for [OII]
; and [SII] corresponding to a density of 100/cm^3 and 10^4 K of:
; <[OII]3726/3729> = 0.75 and <[SII]6716/6731>=1.3 (see IM_TEMDEN for
; details)

; split [OII] 3727 into the [OII] 3726,29 doublet
;   tt = mrdfits(getenv('IMPRO_DIR')+'/etc/temden_table.fits',1)
;   plot, rebin(reform(tt.grid_dens,1,101),185,101), tt.oii_dens_ratio, $
;     ysty=3, psym=8, xr=[1,1E3], /xlog, yr=[0,2]
    doubletratio = 0.75+randomn(seed)*0.1 ; =3726/3729
    factor = [doubletratio/(1.0+doubletratio),1.0/(1.0+doubletratio)]
    this = where(strtrim(flineinfo.name,2) eq '[OII]_3727',comp=therest)
    if this[0] eq -1 then message, 'Missing [OII] 3727!!'
    doublet = im_empty_structure(flineinfo[0],ncopies=2)
    doublet.name = ['[OII]_3726','[OII]_3729']
    doublet.wave = [3726.032D,3728.814D]
    doublet.flux = flineinfo[this].flux*factor
    doublet.ratio = flineinfo[this].ratio*factor
    flineinfo = [flineinfo[therest],doublet]

; add [SII] 6731
    doubletratio = 1.3+randomn(seed)*0.1 ; =6716/6731
    this = where(strtrim(flineinfo.name,2) eq '[SII]_6716')
    if this[0] eq -1 then message, 'Missing [SII] 6716!!'
    doublet = im_empty_structure(flineinfo[0],ncopies=1)
    doublet.name = '[SII]_6731'
    doublet.wave = 6730.815D
    doublet.flux = flineinfo[this].flux/doubletratio
    doublet.ratio = flineinfo[this].ratio/doubletratio
    flineinfo = [flineinfo,doublet]

; add [OIII] 4959
    doubletratio = branch.o_iii ; =5007/4959
    this = where(strtrim(flineinfo.name,2) eq '[OIII]_5007')
    if this[0] eq -1 then message, 'Missing [OIII] 5007!!'
    doublet = im_empty_structure(flineinfo[0],ncopies=1)
    doublet.name = '[OIII]_4959'
    doublet.wave = 4958.911D
    doublet.flux = flineinfo[this].flux/doubletratio
    doublet.ratio = flineinfo[this].ratio/doubletratio
    flineinfo = [flineinfo,doublet]

; add [NII] 6584
    doubletratio = branch.n_ii ; =6584/6548
    this = where(strtrim(flineinfo.name,2) eq '[NII]_6584')
    if this[0] eq -1 then message, 'Missing [NII] 6584!!'
    doublet = im_empty_structure(flineinfo[0],ncopies=1)
    doublet.name = '[NII]_6548'
    doublet.wave = 6548.043D
    doublet.flux = flineinfo[this].flux/doubletratio
    doublet.ratio = flineinfo[this].ratio/doubletratio
    flineinfo = [flineinfo,doublet]
;   struct_print, flineinfo

; combine the hydrogen and forbidden lines
    line = [hlineinfo,flineinfo]
    line = line[sort(line.wave)]
    nline = n_elements(line)
    line.id = lindgen(nline)
;   struct_print, line

; build the wavelength vector and the output spectra; ensure adequate 
; sampling around the lines (use +/-5-sigma)
    if n_elements(wave) eq 0 then begin
       for ii = 0, nline-1 do begin
          width = 4.0*tot_vsigma/light*line[ii].wave
          morewave = range(line[ii].wave-width,line[ii].wave+width,40)
          if ii eq 0 then wave = [line[ii].wave,morewave] else $
            wave = [wave,line[ii].wave,morewave]
       endfor
       wave = wave[sort(wave)]
    endif
    nwave = n_elements(wave)

    if keyword_set(nospectrum) then return, -1
    
    flam_line = wave*0.0
    flam_cont = wave*0.0

; build the emission-line spectrum in log-10 wavelength space, since
; each pixel is a fixed width in velocity space (km/s)
    log10pixsz = 0.2D-4 ; (10^log10pixsz-1)*2.99D5~15 km/s
    minlog10wave = alog10(min(line.wave)*0.9)
    nlog10wave = long((alog10(max(line.wave)*1.1)-minlog10wave)/log10pixsz+1)
    log10wave = minlog10wave + dindgen(nlog10wave)*log10pixsz

    lineflux = log10wave*0.0
    for jj = 0, nline-1 do begin
       amp = line[jj].flux/alog(10)/(sqrt(2.0*!pi)*log10sigma)
       lineflux += amp*exp(-0.5*(log10wave-alog10(line[jj].wave))^2/log10sigma^2)
;      print, line[jj].wave, amp
    endfor
    flam_line = interpolate(lineflux,findex(10.0^log10wave,wave),missing=0.0)
;   djs_plot, 10^log10wave, lineflux, xrange=[4840,4870], ysty=3
;   djs_oplot, wave, flam_line, psym=8, color='red'
    
; now build the nebular continuum spectrum
    flam_cont = nlyc[0]*interpolate(nebcont.flam_neb,findex(nebcont.wave_neb,wave)) ; [erg/s/A]

; scale to a fiducial distance of 10 pc and then add the two
; components together
    factor = 4.0*!dpi*dist^2
    line.flux = line.flux/factor
    
    flam_line = float(flam_line/factor) ; [erg/s/cm^2/A]
    flam_cont = float(flam_cont/factor) ; [erg/s/cm^2/A]
    flam_neb = flam_line + flam_cont
    
return, flam_neb
end
