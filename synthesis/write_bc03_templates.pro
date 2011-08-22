;+
; NAME:
;       WRITE_BC03_TEMPLATES
;
; PURPOSE:
;       Write the BC03 spectra as galaxy templates.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       dlogage - use this input to specify AGEGRID as an argument to
;                 LOGARITHMIC_AGES(); values of 0.1-0.4 are reasonable 
;       agegrid - manually define the ages of the templates [Myr]
;                 (which will be interpolated to the closest matching
;                 SSP age); this keyword is stronger than DLOGAGE
;       suffix  - suffix to append to TFILE, the output template file
;                 name (useful for not overwriting previous templates) 
;       tpath   - output data path for the templates (the default is 
;                 getenv('ISPEC_DIR')+'templates') 
;       minwave - minimum output wavelength (default 1000 A)
;       maxwave - maximum output wavelength (default 25,000 A)
;       dwave   - output wavelength spacing (default 1 A/pixel)
;
; KEYWORD PARAMETERS:
;       salpeter   - read the Salpeter IMF BC03 SSPs (default is the
;                    Chabrier IMF)
;       write      - write out
;
; OUTPUTS:
;       bc   - floating array of all the template [NPIX,NTEMPLATE] 
;       wave - corresponding wavelength vector [NPIX]
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 7, U of A - written
;       jm03sepuofa   - major updates and revisions
;       jm03nov04uofa - more major updates, incorporated IM_READ_BC03() 
;       jm03dec04uofa - added DLOGAGE, AGEGRID, SUFFIX, and TPATH
;                       optional inputs
;       jm04mar03uofa - remove TREMONTI keyword; do not normalize the
;                       templates; added SALPETER keyword; compute M/L
;                       ratios in a large number of bands
;       jm04oct12uofa - use KCORRECT filter routines; use filter
;                       specifications from IM_FILTERSPECS()
;       jm05jul25uofa - output changed to floating-point precision 
;       jm08feb11nyu  - changed default AGEGRID
;       jm08mar13nyu  - added MINWAVE, MAXWAVE, DWAVE optional inputs;
;                       proper calculation of M/L values (now
;                       correctly normalize by the stellar mass at a
;                       given age)
;       jm08aug12nyu  - changed FWHMRES from 3.0 to 3.5 Angstrom
;
; Copyright (C) 2003-2005, 2008, John Moustakas
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

function logarithmic_ages, agemin=agemin, agemax=agemax, $
  dlogage=dlogage, silent=silent
; jm03sep12uofa - return a grid of ages logarithmic in age; all ages
;                 in Myr
; jm08apr19nyu - added AGEMIN, AGEMAX optional inputs    

    if n_elements(agemin) eq 0L then agemin = 1.0D
    if n_elements(agemax) eq 0L then agemax = 13500D
    if n_elements(dlogage) eq 0L then dlogage = 0.20D ; logarithmic age bin

    logage = findgen((alog10(agemax)-alog10(agemin))/dlogage+1.0)*dlogage+alog10(agemin)
    age = 10.0^logage
    age = fix(age[uniq(fix(age))])

    splog, agemin, agemax, dlogage
    
    if (not keyword_set(silent)) then begin
       splog, 'Logarithmic age interval in Myr:'
       niceprint, age
;      niceprint, strjoin(string(age,format='(G0.0)'),', ')
    endif
    
return, age
end    

pro write_bc03_templates, bc, bcwave, bcinfo, dlogage=dlogage, $
  agegrid=agegrid, suffix=suffix, tpath=tpath, minwave=minwave, $
  maxwave=maxwave, dwave=dwave, salpeter=salpeter, write=write

; load filter specifications

    filterlist = [$
      'bessell_'+['U','B','V','R','I'],$
      'sdss_'+['u0','g0','r0','i0','z0'],$
      'twomass_'+['J','H','Ks']]+'.par'
;;  filterlist = ['buser_B3.par','buser_V.par']
    
    filt = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filt)

    sdss = where(filt.sdssflag,nsdss)
    if (nsdss ne 0L) then filt[sdss].vega2ab = 0.0 ; keep the SDSS magnitudes on the AB system

    lsun = 3.826D33 ; [erg/s]
    light = 2.99792458D18       ; speed of light [A/s]
    pc = 3.085678D18            ; conversion from pc to cm [cm]
    fiducial_distance = 10.0*pc ; 10 pc [cm]
    fiducial_mass = 1.0         ; 1 M_sun
    
    if keyword_set(salpeter) then begin
       splog, 'Selecting SALPETER IMF.'
       imf = 'salpeter'
    endif else begin
       splog, 'Selecting CHABRIER IMF.'
       imf = 'chabrier'
    endelse

    logZ = [0.004,0.020,0.050]  ; three metallicity vectors
    strZ = ['Z004','Z02','Z05']
    Zgrid = [2,4,5]
    nZgrid = n_elements(Zgrid)

    if n_elements(tpath) eq 0L then $
      tpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
    if (file_test(tpath,/directory) eq 0L) then begin
       splog, 'Output template path '+tpath+' does not exist.'
       return
    endif

    if n_elements(suffix) eq 0L then suffix = '_' else suffix = '_'+suffix+'_'
    tfile = 'BC03_'+strupcase(strZ)+'_'+imf+suffix+'templates.fits'
    
; choose the SSP ages here

    if n_elements(agegrid) eq 0L then begin

       if n_elements(dlogage) eq 0L then begin
          agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,$
            1500.0,3000.0,7000.0,13000.0]*1E6 & readzero = 0L
;         agegrid = [5,25,100,250,640,1500,12000]*1E6 & readzero = 0L
;;;       agegrid = [0,5,25,100,250,640,1500,12000]*1E6 & readzero = 1L
;;;       agegrid = [0,25,100,250,640,1500,12000]*1E6 & readzero = 1L
;;;       agegrid = [5,25,100,250,640,1500,12000]*1E6 & readzero = 0L

;         agegrid = [5,10,50,100,500,1000,5000,12500]*1E6
;         agegrid = logarithmic_ages(dlogage=0.20,/silent)*1E6
;         agegrid = logarithmic_ages(dlogage=0.35,/silent)*1E6
;         agegrid = [10,25,100,250,640,1500,12000]*1E6
;         agegrid = [1,5,25,100,250,640,1500,12000]*1E6
       endif else begin
          agegrid = logarithmic_ages(dlogage=dlogage,/silent)*1E6
       endelse

    endif else agegrid = agegrid*1E6

    zero = where(agegrid eq 0.0,nzero)
    if (nzero ne 0L) then readzero = 1L else readzero = 0L
    
    nagegrid = n_elements(agegrid)
    
; define the output wavelength range of the templates
    
    fwhmres = 3.0 ; FWHM spectral resolution [Angstrom]
;   fwhmres = 3.5

    if (n_elements(minwave) eq 0L) then minwave = 1000.0  ; minimum wavelength [Angstrom]
    if (n_elements(maxwave) eq 0L) then maxwave = 25000.0 ; maximum wavelength [Angstrom]
    if (n_elements(dwave) eq 0L) then dwave = 1.0         ; [Angstrom/pixel]

    nwave = long((maxwave-minwave)/dwave)+1L
    bcwave = findgen(nwave)*dwave+minwave

; initialize the output info structure

    bcinfo = {$
      templates:                         '', $
      ntemplate:              fix(nagegrid), $
;     template_name:                     '', $
      template_fwhmres:             fwhmres, $
      template_age:                 agegrid, $
      template_Z:                       0.0, $   
      template_imf:                     imf, $
      template_ML_U:       fltarr(nagegrid), $
      template_ML_B:       fltarr(nagegrid), $
      template_ML_V:       fltarr(nagegrid), $
      template_ML_R:       fltarr(nagegrid), $
      template_ML_I:       fltarr(nagegrid)}
;     template_ML_sdss_u:  fltarr(nagegrid), $
;     template_ML_sdss_g:  fltarr(nagegrid), $
;     template_ML_sdss_r:  fltarr(nagegrid), $
;     template_ML_sdss_i:  fltarr(nagegrid), $
;     template_ML_sdss_z:  fltarr(nagegrid), $
;     template_ML_J:       fltarr(nagegrid), $
;     template_ML_H:       fltarr(nagegrid), $
;     template_ML_Ks:      fltarr(nagegrid)}

    for i = 0L, nZgrid-1L do begin

       splog, 'Reading metallicity Z='+string(logZ[i],format='(G0.0)')+'.'

       bcinfo.templates = tfile[i]
       bcgrid = im_read_bc03(metallicity=Zgrid[i],age=agegrid/1E9,salpeter=salpeter,$
         minwave=minwave,maxwave=maxwave,readzero=readzero,/silent,bc03_extras=bce)

       agegrid = bcgrid.age
       bcinfo.template_age = agegrid
;      bcinfo.template_name = 'BC03 ('+strupcase(strZ[i])+')'
       bcinfo.template_Z = logZ[i]

;      bc = fltarr(nwave,nagegrid)
       bc = dblarr(nwave,nagegrid)
       for j = 0L, nagegrid-1L do begin

; COMBINE1FIBER does not work because the dispersion is not constant
          
;         combine1fiber, alog10(bcgrid.wave), bcgrid.flux[*,j], $ ; bcgrid.flux[*,j]*0.0+1.0, $
;           newloglam=alog10(bcwave), newflux=newflux
;         bc[*,j] = lsun*newflux                      ; [erg/s/A/M_sun]

; interpolate the original grid onto our output wavelength vector and
; normalize by the **mass in stars** at the given age

;         stellar_mass = 1.0
          stellar_mass = bce[j].m_
          
          linterp, bcgrid.wave, bcgrid.flux[*,j]/stellar_mass, bcwave, newflux ; [L_sun/A/M_sun]
          bc[*,j] = lsun*newflux                                               ; [erg/s/A/M_sun]
          
; compute the mass-to-light ratio in each filter; Blanton's code needs
; the spectra to be in units of [erg/s/cm2/A], so we create a fiducial
; SED in units of [erg/s/cm2/A] by placing a model having mass
; FIDUCIAL_MASS at a distance FIDUCIAL_DISTANCE 

          knewflux = fiducial_mass*lsun*newflux/(4*!dpi*fiducial_distance^2)   ; [erg/s/cm2/A]

;         bcwave_edges = k_lambda_to_edges(bcwave)
;         fnu = k_project_filters(bcwave_edges,knewflux,filterlist=filterlist) ; [erg/s/cm2/Hz]
          fnu = im_filtermag(bcwave,knewflux,filterlist=filterlist)            ; [erg/s/cm2/Hz]

; compute the apparent and absolute magnitudes in Vega magnitudes and
; in AB magnitudes for the SDSS filters
          
          appmags = -2.5*alog10(fnu) - filt.vega2ab                ; apparent magnitude
          absmags = appmags - (5*alog10(fiducial_distance/pc) - 5) ; absolute=apparent magnitude

          lmratio = 10^(-0.4*(absmags-filt.solarmags)) ; (L/L_sun)/(M/M_sun)
          mlratio = 1.0/lmratio                        ; (M/M_sun)/(L/L_sun)
          
          bcinfo.template_ML_U[j]      = mlratio[0]
          bcinfo.template_ML_B[j]      = mlratio[1]
          bcinfo.template_ML_V[j]      = mlratio[2]
          bcinfo.template_ML_R[j]      = mlratio[3]
          bcinfo.template_ML_I[j]      = mlratio[4]
;         bcinfo.template_ML_sdss_u[j] = mlratio[5]
;         bcinfo.template_ML_sdss_g[j] = mlratio[6]
;         bcinfo.template_ML_sdss_r[j] = mlratio[7]
;         bcinfo.template_ML_sdss_i[j] = mlratio[8]
;         bcinfo.template_ML_sdss_z[j] = mlratio[9]
;         bcinfo.template_ML_J[j]      = mlratio[10]
;         bcinfo.template_ML_H[j]      = mlratio[11]
;         bcinfo.template_ML_Ks[j]     = mlratio[12]

       endfor

;      niceprint, agegrid/1E9, bcinfo.template_ml_b/(bce.m__lb/bce.m_), bcinfo.template_ml_v/(bce.m__lv/bce.m_)
       niceprint, agegrid/1E9, bcinfo.template_ml_b/(bce.m__lb), bcinfo.template_ml_v/(bce.m__lv)

; generate the output header

       mkhdr, header, bc, /extend
       sxdelpar, header, 'COMMENT'
       sxaddpar, header, 'OBJECT', 'BC03 (Z='+string(logZ[i],format='(G0.0)')+') ', ' BC03 model'
       sxaddpar, header, 'CRVAL1', float(minwave), ' central wavelength of first pixel [Angstrom]'
       sxaddpar, header, 'CRPIX1', 1, ' starting pixel (1-indexed)'
       sxaddpar, header, 'CD1_1', float(dwave), ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CTYPE1', 'LINEAR', ' wavelength spacing'
       sxaddpar, header, 'FWHMRES', fwhmres, ' FWHM resolution [Angstrom]'
       sxaddpar, header, 'LOGZ', float(logz[i]), ' metallicity'
       sxaddpar, header, 'IMF', strupcase(imf), ' Initial Mass Function'
       sxaddpar, header, 'DATE', hogg_iso_date(), ' file creation date'

; write out

       if keyword_set(write) and (bcinfo.ntemplate ne 0L) then begin
          splog, 'Writing '+tpath+tfile[i]+'.'
          mwrfits, bcinfo, tpath+tfile[i], /create
          mwrfits, bc, tpath+tfile[i], header
       endif

    endfor 

return
end
