;+
; NAME:
;   BUILD_PEGASE_SSP
;
; PURPOSE:
;   Generate a grid of PEGASE simple stellar populations (SSPs)
;   to be used by BUILD_ISEDFIT_SFHGRID.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   kroupa - use the Kroupa+01 IMF (default is the Salpeter IMF) 
;   cosmic_imf - generate SSPs for my PRIMUS cosmic-IMF project 
;
; OUTPUTS:
;   An information structure is written to getenv('ISEDFIT_SSP_DIR')
;   and the models themselves are written to the 'pegase'
;   subdirectory.
;
; COMMENTS:
;   The code CREATE_PEGASE_SSPS needs to have been run.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Mar 29, UCSD
;
; Copyright (C) 2011, John Moustakas
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

pro build_pegase_ssp, kroupa=kroupa, cosmic_imf=cosmic_imf

    splog, 'Building the PEGASE SSPs'

    pegpath = getenv('PEGASE_HR_DIR')+'/SSPs/'
    ssppath = getenv('ISEDFIT_SSP_DIR')+'/'
    outpath = ssppath+'pegase/'

    imfstr = 'salp'
    if keyword_set(kroupa) then imfstr = 'kroupa01'
    if keyword_set(cosmic_imf) then begin
       alpha2 = primus_cosmicimf_slope() ; high-mass slope
       imfstr = 'cosmicimf_'+string(alpha2,format='(F4.2)')
    endif

    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]
    
; read each SSP in turn, convert to a FITS structure, do some magic,
; and then write out
    Zstr = ['0.0001','0.0004','0.004','0.008','0.02','0.05','0.1']
    nZ = n_elements(Zstr)

    Z = fltarr(nZ)
    sspfile = strarr(nZ)
    for ii = 0, n_elements(imfstr)-1 do begin
       for iZ = 0, nZ-1 do begin
          infile = pegpath+'SSP_'+imfstr[ii]+'_Z'+Zstr[iZ]+'.fits'
          pegase = im_read_peg(infile)
          nage = n_elements(pegase.age)
          npix = n_elements(pegase[0].wave)
          ssp = init_isedfit_ssp(nage=nage,npix=npix)

          ssp.age = pegase.age*1D6 ; [yr]
          ssp.wave = pegase[0].wave
          ssp.flux = pegase.flux
          ssp.mstar = pegase.mstar
          ssp.Z = Zstr[iZ]

; normalize the spectrum by the stellar mass, put it at a fiducial
; distance of 10 pc, and convert to erg/s
;      ssp.flux = ssp.flux/rebin(reform(ssp.mstar,1,nage),npix,nage)
          ssp.flux = ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/Msun]

          sspfile1 = 'pegase_'+imfstr[ii]+'_'+Z2string(ssp.Z)+'.fits'
          im_mwrfits, ssp, outpath+sspfile1, /clobber

          Z[iZ] = ssp.Z
          sspfile[iZ] = sspfile1
       endfor

; write out an information structure
       info = {$
         imf:               imfstr[ii],$
         Z:                      Z,$
         sspfile:          sspfile+'.gz'}

       infofile = ssppath+'info_pegase_'+imfstr[ii]+'.fits'
       im_mwrfits, info, infofile, /clobber
    endfor
    
return
end
