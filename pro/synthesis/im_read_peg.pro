;+
; NAME:
;   IM_READ_PEG()
;
; PURPOSE:
;   Read PEGASE.2 or PEGASE-HR style SSPs into a structure. 
;
; INPUTS:
;   pegfile  - name of the Pegase data file
;
; OPTIONAL INPUTS:
;   mass     - mass normalization
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   peg - output data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;   IM_PEG_STRUCT
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Oct 18, U of A - based on K_READ_PEG,
;     originally written by A. Quintero
;   jm05nov19uofa - added optional mass normalization
;   jm10nov04ucsd - MSTAR should be defined as *all* stars, including
;     remnants, not just *live* stars  
;
; Copyright (C) 2005, 2010, John Moustakas
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

function im_peg_struct, ntime, ncont, nlines
; jm05oct18uofa - compatible with IM_READ_PEG
; jm10nov04ucsd - replaced MALLSTARS with MSTARLIVE
    
    peg = {$
      file:          '',            $
      nage:          ntime,         $
      ncont:         ncont,         $
      nlines:        nlines,        $
      age:           0.0D,          $ ; need double precision
      mgalaxy:       0.0,           $
      mstar:         0.0,           $ ; M_star+M_wd+M_nsbh+M_substellar
      mwd:           0.0,           $ ; white dwarfs
      mnsbh:         0.0,           $ ; neutron stars, black holes
      msubstellar:   0.0,           $
      mstarlive:     0.0,           $ 
      mgas:          0.0,           $
      Zgas:          0.0,           $
      Zstar_mass:    0.0,           $
      Zstar_lbol:    0.0,           $
      lbol:          0.0D,          $ ; double
      tauv:          0.0,           $
      ldust_lbol:    0.0,           $
      sfr:           0.0,           $
      nlyc:          0.0D,          $ ; double
      sniirate:      0.0,           $
      sniarate:      0.0,           $
      starage_mass:  0.0,           $ ; Myr
      starage_lbol:  0.0,           $ ; Myr
      wave:          fltarr(ncont), $
      flux:          dblarr(ncont), $ ; need double precision
      linewave:      fltarr(nlines),$
      lineflux:      dblarr(nlines)}  ; need double precision
    peg = replicate(peg,ntime)

return, peg
end
    
function im_read_peg, pegfile, mass=mass

    if (n_elements(pegfile) eq 0L) then begin
       doc_library, 'im_read_peg'
       return, -1
    endif

    if (file_test(pegfile) eq 0L) then begin
       print, 'Pegase data file '+pegfile+' not found'
       return, -1
    endif
    
    lsun = 3.826D33             ; [erg/s]
    if (n_elements(mass) eq 0L) then mass = 1.0D ; [M_sun]

; read the PEGASE-HR FITS file format; note that if the Basel library
; was used then the wavelength spacing is irregular and the
; wavelengths are stored in the 3rd extension, otherwise the
; wavelength info (ie, from the Elodie library) can be extracted from
; the header (see p. 7 in the Pegase-HR documentation)

    if strmatch(pegfile,'*fits*',/fold) then begin

       hdr = headfits(pegfile)
       ncont = sxpar(hdr,'NAXIS1')
       ntime = sxpar(hdr,'NAXIS2')
       nlines = 61L ; hard-wired!
       peg = im_peg_struct(ntime,ncont,nlines)

       fits_info, pegfile, n_ext=next, /silent

       flux = mrdfits(pegfile,0,head,/silent)
       lines = mrdfits(pegfile,1,/silent)
       info = mrdfits(pegfile,2,/silent)

       if (next eq 3) then begin
          waveinfo = mrdfits(pegfile,3,/silent)
          wave = waveinfo.bfit
       endif else wave = make_wave(hdr)
       
       peg.file     = file_basename(pegfile)
       peg.wave     = wave
       peg.flux     = flux*mass*lsun
       peg.linewave = lines.wave
       peg.lineflux = transpose(lines.fluxline*mass*lsun)

       peg.age         = info.age
       peg.mgalaxy     = info.mgal*mass
       peg.mwd         = info.mwd*mass
       peg.mnsbh       = info.mbhns*mass
       peg.msubstellar = info.msub*mass
       peg.mstarlive   = info.mstars*mass ; all alive stars
       peg.mstar       = peg.mstarlive+peg.mwd+peg.mnsbh+peg.msubstellar ; all stars

       peg.mgas        = info.sigmagas*mass
       peg.Zgas        = info.zgas; *info.sigmagas ; gas-phase metallicity
       peg.Zstar_mass  = info.zstars
       peg.Zstar_lbol  = info.zbol

       peg.lbol         = info.fluxbol*mass*lsun
       peg.tauv         = info.tauv
       peg.ldust_lbol   = info.fluxext
       peg.sfr          = info.sfr*mass/1D6 ; [M_sun/yr]
       peg.nlyc         = info.nlymtot*mass
       peg.sniirate     = info.nsniitot*mass
       peg.sniarate     = info.nsniatot*mass
       peg.starage_mass = info.agestars
       peg.starage_lbol = info.agebol

       return, peg
    endif
    
    openr, unit, pegfile, /get_lun

    a = ' '
    readf, unit, a
    while a ne '************************************************************' $
      do begin
       readf, unit, a
    endwhile

    point = lonarr(3)
    readf, unit, point 

    lam = lonarr(point[1])
    readf, unit, lam

    lines = lonarr(point[2])
    readf, unit, lines

    arr1 = dblarr(10) ; note double arrays!
    arr2 = dblarr(9)
    cont_lum_arr = dblarr(point[1])
    line_lum_arr = dblarr(point[2])

    peg = im_peg_struct(point[0], point[1], point[2])

    for i = 0L, point[0]-1L do begin

       readf, unit, arr1 
       readf, unit, arr2
       readf, unit, cont_lum_arr
       readf, unit, line_lum_arr

       peg.file        = file_basename(pegfile)

       peg.wave        = lam
       peg[i].flux     = cont_lum_arr*mass
       peg.linewave    = lines
       peg[i].lineflux = line_lum_arr*mass*lsun

       peg[i].age         = arr1[0]
       peg[i].mgalaxy     = arr1[1]*mass
       peg[i].mwd         = arr1[3]*mass
       peg[i].mnsbh       = arr1[4]*mass
       peg[i].msubstellar = arr1[5]*mass 
       peg[i].mstarlive   = arr1[2]*mass
       peg[i].mstar       = peg[i].mstarlive+peg[i].mwd+peg[i].mnsbh+peg[i].msubstellar

       peg[i].mgas        = arr1[6]*mass
       peg[i].Zgas        = arr1[7]
       peg[i].Zstar_mass  = arr1[8]
       peg[i].Zstar_lbol  = arr1[9]

       peg[i].lbol         = arr2[0]*mass
       peg[i].tauv         = arr2[1]
       peg[i].ldust_lbol   = arr2[2]
       peg[i].sfr          = arr2[3]*mass/1D6 ; [M_sun/yr]
       peg[i].nlyc         = arr2[4]*mass
       peg[i].sniirate     = arr2[5]*mass
       peg[i].sniarate     = arr2[6]*mass
       peg[i].starage_mass = arr2[7]
       peg[i].starage_lbol = arr2[8]

    endfor

    free_lun, unit
    
return, peg    
end
