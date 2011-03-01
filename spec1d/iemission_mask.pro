;+
; NAME:
;   IEMISSION_MASK()
;
; PURPOSE:
;   Mask various UV/optical spectral features, including nebular and
;   QSO emission lines, sky lines, and the telluric bandpasses.
;
; INPUTS:
;   wave - *observed* frame wavelength vector [Angstrom] 
;
; OPTIONAL INPUTS:
;   z     - redshift (default 0.0)
;   vdisp - velocity dispersion (default 150 km/s)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   mask  - pixel mask for WAVE (1=good, 0=masked)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 12, NYU - written, based on EMISSION_MASK() 
;
; Copyright (C) 2008, John Moustakas
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

function iemission_mask, wave, z=z, vdisp=vdisp, nebular=nebular, $
  qso=qso, sky=sky, telluric=telluric

    nwave = n_elements(wave)
    if (nwave eq 0L) then begin
       doc_library, 'emission_mask'
       return, -1
    endif

    light = 2.99792458D5 ; speed of light [km/s]
    if (n_elements(z) eq 0L) then z = 0.0
    if (n_elements(vdisp) eq 0L) then vdisp = 150.0 ; [km/s]
    if (n_elements(nebfactor) eq 0L) then nebfactor = 150.0 ; [km/s]
    
    restwave = wave/(1.0+z)
    
; mask the nebular lines    
    
    neblinelist = [$
      3425.868,$                ; [OII]
      3703.860,$                ; He I
      3726.032,$                ; [OII]
      3728.815,$                ; [OII]
      3750.151,$                ; H12
      3770.630,$                ; H11
      3797.898,$                ; H10
      3819.64, $                ; He I
      3835.384,$                ; H9
      3868.75, $                ; [NeIII]
      3889.049,$                ; H8
      3970.072,$                ; H-epsilon
      4026.21, $                ; He I
      4068.60, $                ; [S II]
      4101.734,$                ; H-delta
      4340.46 ,$                ; H-gamma
      4363.21 ,$                ; [OIII]
      4471.50, $                ; He I
      4861.33 ,$                ; H-beta
      4958.91 ,$                ; [OIII]
      5006.84 ,$                ; [OIII]
      5200.26, $                ; [N I]
      5875.96 ,$                ; HeI
      5890.0  ,$                ; Na D doublet
      5896.0  ,$                ; Na D doublet
      6300.30 ,$                ; OI
      6312.40, $                ; [S III]
      6363.78, $                ; [O I]
      6548.04 ,$                ; NII
      6562.80 ,$                ; H-alpha
      6583.46 ,$                ; NII
      6678.15, $                ; He I
      6716.14, $                ; SII
      6730.81, $                ; SII
      7065.28, $                ; He I
      7135.78, $                ; [Ar III]
      7319.65, $                ; [O II]
      7330.16, $                ; [O II]
      7751.12  $                ; [Ar III]
      ]

    nebfactor = 4.0
    nebmask_width = nebfactor*(100.0>vdisp<500.0)
    nebmask = make_array(nwave,/byte,value=1)

    if keyword_set(nebular) then begin
       for ii = 0L, n_elements(neblinelist)-1L do $
         nebmask = nebmask and (abs(light*(restwave-neblinelist[ii])/$
         neblinelist[ii]) gt nebmask_width)
    endif
    
; mask QSO emission lines
    
    qsolinelist = [$
      1218.0,$
      1305.0,$
      1400.0,$
      1549.0,$  ; [CIV]
      1909.0,$  ; CIII]
      2799.495$ ; Mg II
;     4861.33,$                ; H-beta
;     6562.80$                ; H-alpha
    ]

    qsofactor = 10.0
    qsomask_width = qsofactor*(1000.0 > vdisp < 5000.0)
    qsomask = make_array(nwave,/byte,value=1)

    if keyword_set(qso) then begin
       for ii = 0L, n_elements(qsolinelist)-1L do $
         qsomask = qsomask and (abs(light*(restwave-qsolinelist[ii])/$
         qsolinelist[ii]) gt qsomask_width)
    endif

; mask strong sky lines; use NEBMASK_WIDTH and note that we use the
; *observed* wavelengths here
    
    skylinelist = [5577.339,5889.950,6300.32,6363.81] ; sky wavelengths

    skymask = make_array(nwave,/byte,value=1)
    skymask_width = nebmask_width ; note!

    if keyword_set(sky) then begin
       for ii = 0L, n_elements(skylinelist)-1L do $
         skymask = skymask and (abs(light*(wave-skylinelist[ii])/$
         skylinelist[ii]) gt skymask_width)
    endif

; finally flag telluric features and the add all the masks together 

    tellmask = wave*0.0+1.0
    if keyword_set(telluric) then begin
       tellmask = telluric_mask(wave)
    endif

    outmask = (nebmask and qsomask) and (skymask and tellmask)

return, outmask
end
