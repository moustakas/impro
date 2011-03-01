;+
; NAME:
;       WRITE_ANCILLARY_DATA
;
; PURPOSE:
;       Compile and write ancillary data for an input set of galaxies.  
;
; INPUTS:
;       datapath  - 
;       outpath   - 
;       basicname - 
;       photoname - 
;       distname  - 
;       diamname  - 
;       leda      - 
;       outname   - 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       norasort - do not sort by RA
;       write    - write out
;
; OUTPUTS:
;       table - output data table
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 May 16 - generalized from earlier code 
;       jm05jun02uofa - use inclination angles computed from R25
;                       assuming a fixed intrinsic axial ratio
;                       (cf. Tully et al. 1998) rather than the LEDA
;                       inclination angles; compute the extinction
;                       correction to face-on inclination
;       jm05jul21uofa - replaced the LEDANAME optional input with
;                       LEDA, which is now an optional input data 
;                       structure 
;       jm05jul22uofa - removed FIXRC3ATLAS keyword and RC3RAD
;                       optional input; added DIAMNAME optional input
;       jm05aug17uofa - require an RC3 B magnitude when collecting UV
;                       RC3 magnitudes (see comments, below)
;       jm06apr04uofa - added NORASORT
;       jm07nov29nyu  - added Sanders et al. IRAS fluxes
;       jm08jun06nyu  - a bit of a major rewrite; consolidate the
;                       photometry into single K-correct-like arrays
;                       (see also INIT_ANCILLARY_DATA); remove some of
;                       the LEDA stuff, and also some of the old radio
;                       and UV fluxes and other ancillary measurements
;                       that I have *never* used!
;
; Copyright (C) 2005-2008, John Moustakas
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

pro write_ancillary_data, table, datapath=datapath, outpath=outpath, $
  basicname=basicname, photoname=photoname, distname=distname, $
  diamname=diamname, leda=leda, outname=outname, norasort=norasort, $
  write=write

; initialize path names    

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = cwd()

    nbasic = n_elements(basicname)
    nphoto = n_elements(photoname)
    ndist = n_elements(distname)
    ndiam = n_elements(diamname)
    nleda = n_elements(leda)
    
    if (nbasic eq 0L) then begin
       splog, 'BASICNAME must be provided.'
       return
    endif

;   if (nbasic eq 0L) or (nphoto eq 0L) or (ndist eq 0L) or (ndiam eq 0L) then begin
;      splog, 'BASICNAME, PHOTONAME, DISTNAME, and DIAMNAME must be provided.'
;      return
;   endif

    if keyword_set(write) and (n_elements(outname) eq 0L) then begin
       splog, 'OUTNAME must be provided.'
       return
    endif

; define the cosmology and some constants    
    
    red, omega0=0.3, omegalambda=0.7, h100=0.70
    h100 = redh100()
    
    mpc2m = 3.086D22     ; [m/Mpc]
    mpc2cm = 3.086D24    ; [cm/Mpc]
    ergs2watts = 1D7     ; [(erg/s)/W]
    light = 2.99792458D5 ; speed of light [km/s]
    lsun = 3.826D33      ; bolometric solar luminosity [erg/s]
    mbolsun = 4.74       ; bolometric absolute solar magnitude [mag]
    
; compile some filter information    
    
    filterlist = [$
      'bessell_U',$
      'bessell_B',$
      'bessell_V',$
      'bessell_R',$
      'bessell_I',$
      'twomass_J',$
      'twomass_H',$
      'twomass_Ks'$
      ]+'.par'

    filtinfo = im_filterspecs(filterlist=filterlist)

    Bmsun = filtinfo[1].solarmags

    U_weff = filtinfo[0].weff
    B_weff = filtinfo[1].weff
    V_weff = filtinfo[2].weff
    R_weff = filtinfo[3].weff
    I_weff = filtinfo[4].weff

    J_weff = filtinfo[5].weff
    H_weff = filtinfo[6].weff
    K_weff = filtinfo[7].weff
    
; retrieve the NED basic data, NED photometry, distances, and NED
; diameters 

    if (file_test(datapath+basicname) eq 0L) then begin
       splog, 'File '+datapath+basicname+' not found.'
       return
    endif 

    basicdata = mrdfits(datapath+basicname,1,/silent)    
    ngalaxy = n_elements(basicdata)

    if (nleda ne 0L) then begin
       if (nleda ne ngalaxy) then begin
          splog, 'BASICDATA and LEDA must be the same size.'
          return
       endif
;      niceprint, basicdata.galaxy, photodata.galaxy, leda.name
    endif; else niceprint, basicdata.galaxy, photodata.galaxy

; read other ancillary data    
    
    if (nphoto ne 0L) then begin
       if (file_test(datapath+photoname) eq 0L) then begin
          splog, 'File '+datapath+photoname+' not found.'
          return
       endif else begin
          photodata = mrdfits(datapath+photoname,1,/silent)
          if (ngalaxy ne n_elements(photodata)) then begin
             splog, 'BASICDATA and PHOTODATA must be the same size.'
             return
          endif
       endelse
    endif

    if (ndiam ne 0L) then begin
       if (file_test(datapath+diamname) eq 0L) then begin
          splog, 'File '+datapath+diamname+' not found.'
          return
       endif else begin
          diamdata = mrdfits(datapath+diamname,1,/silent)
          if (ngalaxy ne n_elements(diamdata)) then begin
             splog, 'BASICDATA and DIAMDATA must be the same size.'
             return
          endif
       endelse
    endif

    if (ndist ne 0L) then begin
       if (file_test(datapath+distname) eq 0L) then begin
          splog, 'File '+datapath+distname+' not found.'
          return
       endif else begin
          distdata = mrdfits(datapath+distname,1,/silent)
          if (ngalaxy ne n_elements(distdata)) then begin
             splog, 'BASICDATA and DISTDATA must be the same size.'
             return
          endif
       endelse 
    endif
    
; initialize the output data table

    table = init_ancillary_data(ngalaxy=ngalaxy)
    
; fill the structure with preliminaries

    if (nleda ne 0L) then table.leda_galaxy = 'PGC'+string(leda.pgc,format='(I7.7)')
    
    table.galaxy         = strupcase(strcompress(basicdata.galaxy,/remove))
    table.ned_galaxy     = strupcase(strcompress(basicdata.ned_galaxy,/remove))
    table.ra             = basicdata.ra
    table.dec            = basicdata.dec
    table.z              = basicdata.z
    table.z_err          = basicdata.z_err
    table.ned_morphology = strcompress(basicdata.morph,/remove)

    zgood = where(table.z gt -900.0,nzgood)
    if (nzgood ne 0L) then table[zgood].cz = table[zgood].z*light
    zgood = where(table.z_err gt -900.0,nzgood)
    if (nzgood ne 0L) then table[zgood].cz_err = table[zgood].z_err*light

    if (ndist ne 0L) then begin
       table.distance        = distdata.distance
       table.distance_err    = distdata.distance_err
       table.distance_ref    = distdata.distance_ref
       table.distance_texref = distdata.distance_texref
       table.distance_method = distdata.distance_method
       table.distance_model  = distdata.modeldist
    endif

    ra = 15.0D*im_hms2dec(table.ra)
    dec = im_hms2dec(table.dec)

; NED classification and special notes

    table.ned_class = strtrim(basicdata.class,2)
    table.ned_note = strtrim(basicdata.note,2)
    
; determine the E(B-V) color excess in our Galaxy along the line of
; sight using the SFD dust maps.  convert RA and DEC to Galactic
; coordinates.

    splog, 'Computing the foreground Galactic extinction values.'

    glactc, 15.0D*im_hms2dec(table.ra), im_hms2dec(table.dec), 2000.0, gl, gb, 1, /degree
    table.ebv_mw = dust_getval(gl,gb,/interp)

; compute the physical scale [pc/arcsec]

    good = where((table.distance gt -900.0),ngood)
    if (ngood ne 0L) then begin
       dz = h100*100.0*table[good].distance/light      ; cosmological redshift
       table[good].xscale = dangular(dz,/pc)/206265.0D ; [pc/arcsec]
;      niceprint, table[good].galaxy, dz, table[good].distance, table[good].xscale
    endif
    
; RC3 magnitudes from NED; do not use "corrected" (e.g., B_T^0) RC3
; magnitudes, only "uncorrected" magnitudes (e.g., B_T or m_B);
; UPDATE: generate a mini-table here with all the quantities of
; interest, and then concatenate them into TABLE below (jm08jun06nyu) 

    if (nphoto ne 0L) then begin

       splog, 'Parsing optical photometry from NED/LEDA.'

       minitable = {$
         rc3_u_ref:                '...', $ ; RC3 or LEDA?
         rc3_u:                      -999.0, $ ; RC3/LEDA total apparent U magnitude and error
         rc3_u_err:                  -999.0, $
         rc3_b_ref:                '...', $ ; RC3 or LEDA?
         rc3_b:                      -999.0, $ ; RC3/LEDA total apparent B magnitude and error
         rc3_b_err:                  -999.0, $
         rc3_v_ref:                '...', $ ; RC3 or LEDA?
         rc3_v:                      -999.0, $ ; RC3/LEDA total apparent V magnitude and error
         rc3_v_err:                  -999.0, $
         ebv_mw:                     -999.0}
       minitable = replicate(minitable,ngalaxy)
       minitable.ebv_mw = table.ebv_mw
       
       good = where((strmatch(strtrim(photodata.B_RC3_FLAG,2),'B (B_T)',/fold) eq 1B) or $
         ((strmatch(strtrim(photodata.B_RC3_FLAG,2),'B (m_B)') eq 1B)),ngood)
       if ngood ne 0L then begin
          minitable[good].RC3_B_ref = 'RC3/'+strtrim(photodata[good].b_rc3_flag,2)
          minitable[good].RC3_B = photodata[good].B_RC3
          minitable[good].RC3_B_err = photodata[good].B_RC3_err
       endif
;      niceprint, photodata[good].b_rc3, photodata[good].b_rc3_err, photodata[good].b_rc3_flag

; only use U- and V-band magnitudes if the B-band magnitude is also
;             defined; I have found (05aug17) that occasionally UV
;             magnitudes without a B magnitude are unreliable (e.g,
;             ARP256N=ARP256NED02)
       
       good = where((strmatch(strtrim(photodata.U_rc3_FLAG,2),'U (U_T)',/fold) eq 1B) and $
         (minitable.rc3_b gt -900.0),ngood)
       if ngood ne 0L then begin
          minitable[good].RC3_U_ref = 'RC3/'+strtrim(photodata[good].u_rc3_flag,2)
          minitable[good].RC3_U = photodata[good].U_RC3
          minitable[good].RC3_U_err = photodata[good].U_RC3_err
       endif
;      niceprint, photodata[good].u_rc3, photodata[good].u_rc3_err, photodata[good].u_rc3_flag

       good = where((strmatch(strtrim(photodata.V_RC3_FLAG,2),'V (V_T)',/fold) eq 1B) and $
         (minitable.rc3_b gt -900.0),ngood)
       if ngood ne 0L then begin
          minitable[good].RC3_V_ref = 'RC3/'+strtrim(photodata[good].v_rc3_flag,2)
          minitable[good].RC3_V = photodata[good].V_RC3
          minitable[good].RC3_V_err = photodata[good].V_RC3_err
       endif
;      niceprint, photodata[good].v_rc3, photodata[good].v_rc3_err, photodata[good].v_rc3_flag

       if (nleda ne 0L) then begin

          if tag_exist(leda,'bt') then begin
             need = where((minitable.rc3_b eq -999.0) and (strcompress(leda.bt,/remove) ne ''),nneed)
             if (nneed ne 0L) then begin
                minitable[need].rc3_B_ref = 'LEDA'
                minitable[need].rc3_B = leda[need].bt
                minitable[need].rc3_B_err = leda[need].e_bt
                splog, 'Gathering LEDA B photometry for '+string(nneed,format='(I0)')+' galaxies.'
             endif
          endif

          if tag_exist(leda,'bt') and tag_exist(leda,'ubt') then begin
             need = where((minitable.rc3_U eq -999.0) and (strcompress(leda.bt,/remove) ne '') and $
               (strcompress(leda.ubt,/remove) ne '') and (strcompress(leda.ubt,/remove) ne '0'),nneed)
             if (nneed ne 0L) then begin
                minitable[need].rc3_U_ref = 'LEDA'
                minitable[need].rc3_U = float(leda[need].ubt)+float(leda[need].bt)
                minitable[need].rc3_U_err = leda[need].e_bt
                splog, 'Gathering LEDA U photometry for '+string(nneed,format='(I0)')+' galaxies.'
             endif
          endif

          if tag_exist(leda,'bt') and tag_exist(leda,'bvt') then begin
             need = where((minitable.rc3_V eq -999.0) and (strcompress(leda.bt,/remove) ne '') and $
               (strcompress(leda.bvt,/remove) ne '') and (strcompress(leda.bvt,/remove) ne '0'),nneed)
             if (nneed ne 0L) then begin
                minitable[need].rc3_V_ref = 'LEDA'
                minitable[need].rc3_V = float(leda[need].bt)-float(leda[need].bvt)
                minitable[need].rc3_V_err = leda[need].e_bt
                splog, 'Gathering LEDA V photometry for '+string(nneed,format='(I0)')+' galaxies.'
             endif
          endif

;;; mean B-band surface brightness
;;          
;;          if tag_exist(leda,'bri25') then begin
;;             need = where(strcompress(leda.bri25,/remove) ne '',nneed)
;;             if (nneed ne 0L) then table[need].sb25_b = leda[need].bri25
;;          endif
          
       endif

;      need = where(table.rc3_u lt -900.0,nneed)
;      if (nneed ne 0L) then begin
;         splog, 'The following objects need U-band photometry.'
;         niceprint, table[need].galaxy
;      endif
          
;      need = where(table.rc3_b lt -900.0,nneed)
;      if (nneed ne 0L) then begin
;         splog, 'The following objects need B-band photometry.'
;         niceprint, table[need].galaxy
;      endif
       
;      need = where(table.rc3_v lt -900.0,nneed)
;      if (nneed ne 0L) then begin
;         splog, 'The following objects need V-band photometry.'
;         niceprint, table[need].galaxy
;      endif

; correct the UV/optical photometric data for Galactic extinction

       good = where(minitable.rc3_U gt -900,ngood)
       if ngood ne 0L then minitable[good].rc3_U = minitable[good].rc3_U - $
         minitable[good].ebv_mw*k_lambda(U_weff,/odonnell)
       good = where(minitable.rc3_B gt -900,ngood)
       if ngood ne 0L then minitable[good].rc3_B = minitable[good].rc3_B - $
         minitable[good].ebv_mw*k_lambda(B_weff,/odonnell)
       good = where(minitable.rc3_V gt -900,ngood)
       if ngood ne 0L then minitable[good].rc3_V = minitable[good].rc3_V - $
         minitable[good].ebv_mw*k_lambda(V_weff,/odonnell)

; now merge MINITABLE into TABLE

       table.rc3_ubv_ref = transpose([[minitable.rc3_u_ref],$
         [minitable.rc3_b_ref],[minitable.rc3_v_ref]])
       table.rc3_ubv = transpose([[minitable.rc3_u],$
         [minitable.rc3_b],[minitable.rc3_v]])
       table.rc3_ubv_err = transpose([[minitable.rc3_u_err],$
         [minitable.rc3_b_err],[minitable.rc3_v_err]])

       for iband = 0L, 2L do begin
          good = where((table.distance gt -900.0) and (table.rc3_ubv[iband] gt -900.0),ngood)
          if (ngood ne 0L) then begin
             dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
               err_dist=table[good].distance_err)
             table[good].rc3_ubv_absmag[iband] = table[good].rc3_ubv[iband] - dmod
             table[good].rc3_ubv_absmag_err[iband] = sqrt(table[good].rc3_ubv_err[iband]^2.0 + err_dmod^2.0)
          endif
       endfor

    endif 
       
; 2MASS near-infrared photometry.  correct for Galactic extinction.
; we use the Cardelli, Clayton, & Mathis (1989) extinction curve and
; the O'Donnell (1994) coefficients; ; UPDATE: generate a mini-table
; here with all the quantities of interest, and then concatenate them
; into TABLE below (jm08jun06nyu)


    if (nphoto ne 0L) then begin
       
       splog, 'Parsing 2MASS photometry.'

; total magnitudes
       
       minitable = {$
         twomass_JHK_ref:              '...', $ ; XSC or LGA?
         twomass_J:                      -999.0, $ ; 2MASS J magnitude and error
         twomass_J_err:                  -999.0, $
         twomass_H:                      -999.0, $ ; 2MASS H magnitude and error
         twomass_H_err:                  -999.0, $
         twomass_K:                      -999.0, $ ; 2MASS Ks magnitude and error
         twomass_K_err:                  -999.0, $
         ebv_mw:                         -999.0}
       minitable = replicate(minitable,ngalaxy)
       minitable.ebv_mw = table.ebv_mw
       
       minitable.twomass_jhk_ref = photodata.J_2MASS_ref
       
       good = where(photodata.J_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_J  = photodata[good].J_2mass - $
         minitable[good].ebv_mw*k_lambda(J_weff,/odonnell)
       good = where(photodata.J_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_J_err  = photodata[good].J_2mass_err

       good = where(photodata.H_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_H  = photodata[good].H_2mass - $
         minitable[good].ebv_mw*k_lambda(H_weff,/odonnell)
       good = where(photodata.H_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_H_err  = photodata[good].H_2mass_err

       good = where(photodata.K_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_K  = photodata[good].K_2mass - $
         minitable[good].ebv_mw*k_lambda(K_weff,/odonnell)
       good = where(photodata.K_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_K_err  = photodata[good].K_2mass_err

; merge MINITABLE into TABLE and calculate absolute magnitudes

       table.twomass_jhk_ref = minitable.twomass_jhk_ref
       table.twomass_jhk = transpose([[minitable.twomass_j],$
         [minitable.twomass_h],[minitable.twomass_k]])
       table.twomass_jhk_err = transpose([[minitable.twomass_j_err],$
         [minitable.twomass_h_err],[minitable.twomass_k_err]])

       for iband = 0L, 2L do begin
          good = where((table.distance gt -900.0) and (table.twomass_jhk[iband] gt -900.0),ngood)
          if (ngood ne 0L) then begin
             dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
               err_dist=table[good].distance_err)
             table[good].twomass_jhk_absmag[iband] = table[good].twomass_jhk[iband] - dmod
             table[good].twomass_jhk_absmag_err[iband] = sqrt(table[good].twomass_jhk_err[iband]^2.0 + err_dmod^2.0)
          endif
       endfor

; isophotal magnitudes

       minitable = {$
         twomass_JHK20_ref:              '...', $ ; XSC or LGA?
         twomass_J20:                      -999.0, $ ; 2MASS J magnitude and error
         twomass_J20_err:                  -999.0, $
         twomass_H20:                      -999.0, $ ; 2MASS H magnitude and error
         twomass_H20_err:                  -999.0, $
         twomass_K20:                      -999.0, $ ; 2MASS Ks magnitude and error
         twomass_K20_err:                  -999.0, $
         ebv_mw:                         -999.0}
       minitable = replicate(minitable,ngalaxy)
       minitable.ebv_mw = table.ebv_mw
       
       minitable.twomass_jhk20_ref = photodata.J20_2MASS_ref
       
       good = where(photodata.J20_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_J20  = photodata[good].J20_2mass - $
         minitable[good].ebv_mw*k_lambda(J_weff,/odonnell)
       good = where(photodata.J20_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_J20_err  = photodata[good].J20_2mass_err

       good = where(photodata.H20_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_H20  = photodata[good].H20_2mass - $
         minitable[good].ebv_mw*k_lambda(H_weff,/odonnell)
       good = where(photodata.H20_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_H20_err  = photodata[good].H20_2mass_err

       good = where(photodata.K20_2mass gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_K20  = photodata[good].K20_2mass - $
         minitable[good].ebv_mw*k_lambda(K_weff,/odonnell)
       good = where(photodata.K20_2mass_err gt -900.0,ngood)
       if ngood ne 0L then minitable[good].twomass_K20_err  = photodata[good].K20_2mass_err

; merge MINITABLE into TABLE and calculate absolute magnitudes

       table.twomass_jhk20 = transpose([[minitable.twomass_j20],$
         [minitable.twomass_h20],[minitable.twomass_k20]])
       table.twomass_jhk20_err = transpose([[minitable.twomass_j20_err],$
         [minitable.twomass_h20_err],[minitable.twomass_k20_err]])

       for iband = 0L, 2L do begin
          good = where((table.distance gt -900.0) and (table.twomass_jhk20[iband] gt -900.0),ngood)
          if (ngood ne 0L) then begin
             dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
               err_dist=table[good].distance_err)
             table[good].twomass_jhk20_absmag[iband] = table[good].twomass_jhk20[iband] - dmod
             table[good].twomass_jhk20_absmag_err[iband] = sqrt(table[good].twomass_jhk20_err[iband]^2.0 + err_dmod^2.0)
          endif
       endfor

    endif 

; IRAS fluxes and errors: prioritize according to the following order
; (see, e.g., Bell 2003): Sanders et al. 2003; Rice et al. 1988, ApJS,
; 68, 91 for large optical galaxies; Soifer et al. 1989, AJ, 98, 766
; for the Bright Galaxy Sample; and Moshir et al. 1990 for the IRAS
; Faint Source Catalog; MAJOR UPDATE: this prioritization is now done
; within NED_WEBGET_PHOTO! (jm08jan15nyu)

    if (nphoto ne 0L) then begin

       splog, 'Parsing IRAS photometry.'

       table.iras_ref = transpose([[photodata.iras_12_ref],$
         [photodata.iras_25_ref],[photodata.iras_60_ref],[photodata.iras_100_ref]])
       table.iras = transpose([[photodata.iras_12],[photodata.iras_25],$
         [photodata.iras_60],[photodata.iras_100]])
       table.iras_err = transpose([[photodata.iras_12_err],$
         [photodata.iras_25_err],[photodata.iras_60_err],[photodata.iras_100_err]])

; for objects without good 12 or 25 micron fluxes (including upper
; limits), use the relations in Bell (2003) to predict them, using the
; 100- and 60-micron fluxes, respectively; the error in these
; predictions is ~30%

       predict_error = 0.30     ; [%]
       
       good = where(table.iras[3] gt 0.0,ngood)
       if (ngood ne 0L) then begin
          table[good].iras_12_predict = 0.0326*table[good].iras[3]
          table[good].iras_12_predict_err = table[good].iras_12_predict*predict_error
       endif
       
       good = where(table.iras[2] gt 0.0,ngood)
       if (ngood ne 0L) then begin
          table[good].iras_25_predict = 0.131*table[good].iras[2]
          table[good].iras_25_predict_err = table[good].iras_12_predict*predict_error
       endif

; now compute absolute magnitudes       
       
       for iband = 0L, 3L do begin
          good = where((table.distance gt -900.0) and (table.iras[iband] gt -900.0),ngood)
          if (ngood ne 0L) then begin
             dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
               err_dist=table[good].distance_err)
             table[good].iras_absmag[iband] = table[good].iras[iband] - dmod
             table[good].iras_absmag_err[iband] = sqrt(table[good].iras_err[iband]^2.0 + err_dmod^2.0)
          endif
       endfor

    endif 

; morphological types from LEDA

    if (nleda ne 0L) then begin
    
       splog, 'Parsing morphological types from LEDA.'

; morphological type

       if tag_exist(leda,'type') then begin
          need = where((strcompress(table.rc3_type,/remove) eq '') and (strcompress(leda.type,/remove) ne ''),nneed)
          if (nneed ne 0L) then table[need].rc3_type = leda[need].type
       endif
       
       if tag_exist(leda,'t') then begin
          need = where((table.rc3_t eq -999.0) and (strcompress(leda.t,/remove) ne ''),nneed)
          if (nneed ne 0L) then table[need].rc3_t = leda[need].t
       endif

; visual morphological characteristics
       
;      if tag_exist(leda,'bar') then begin
;         bar = where(strcompress(leda.bar,/remove) eq 'B',nbar)
;         if nbar ne 0L then table[bar].bar = 1L
;      endif
;
;      if tag_exist(leda,'ring') then begin
;         ring = where(strcompress(leda.ring,/remove) eq 'R',nring)
;         if nring ne 0L then table[ring].ring = 1L
;      endif
;
;      if tag_exist(leda,'multiple') then begin
;         multiple = where(strcompress(leda.multiple,/remove) eq 'M',nmultiple)
;         if nmultiple ne 0L then table[multiple].multiple = 1L
;      endif

    endif
    
; diameters

    if (ndiam ne 0L) then begin

       splog, 'Parsing diameters and position angles from NED/LEDA.'

; 2MASS diameters    

       good = where((diamdata.twomass_major_axis gt -900.0) and $
         (diamdata.twomass_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].twomass_maj = diamdata[good].twomass_major_axis/60.0
          table[good].twomass_min = diamdata[good].twomass_minor_axis/60.0
          table[good].twomass_diameter_ref = diamdata[good].twomass_reference
       endif

       good = where((diamdata.twomass_k20_major_axis gt -900.0) and $
         (diamdata.twomass_k20_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].twomass_k20_maj = diamdata[good].twomass_k20_major_axis/60.0
          table[good].twomass_k20_min = diamdata[good].twomass_k20_minor_axis/60.0
;         table[good].twomass_diameter_ref = diamdata[good].twomass_reference ; assumed the same as above
       endif

       good = where((diamdata.twomass_posangle gt -900.0),ngood)
       if (ngood ne 0L) then table[good].twomass_posangle = diamdata[good].twomass_posangle

; RC3 diameters    
       
       good = where((table.d25_maj lt -900.0) and (diamdata.rc3_major_axis gt -900.0) and $
         (diamdata.rc3_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].d25_maj = diamdata[good].rc3_major_axis/60.0
          table[good].d25_min = diamdata[good].rc3_minor_axis/60.0
          table[good].d25_ref = 'RC3'
          table[good].d25_ref = diamdata[good].rc3_reference
;         struct_print, struct_trimtags(diamdata[good],select=['GALAXY','RC3*'])
       endif

       good = where((table.optical_posangle lt -900.0) and (diamdata.rc3_posangle gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].optical_posangle = diamdata[good].rc3_posangle
          table[good].optical_posangle_ref = 'RC3'
          table[good].optical_posangle_ref = diamdata[good].rc3_reference
       endif

; LEDA

       if (nleda ne 0L) then begin

          if tag_exist(leda,'logd25') then begin
             need = where((table.d25_maj lt -900.0) and (strcompress(leda.logd25,/remove) ne ''),nneed)
             if (nneed ne 0L) then begin
                table[need].d25_maj = 0.1*10.0^float(leda[need].logd25)
                table[need].d25_ref = 'LEDA'
                table[need].d25_ref = 'LEDA'
             endif
          endif

          if tag_exist(leda,'logr25') then begin
             need = where((table.d25_min lt -900.0) and (table.d25_maj gt -900.0) and $
               (strcompress(leda.logr25,/remove) ne ''),nneed)
             if (nneed ne 0L) then table[need].d25_min = 10.0^(-float(leda[need].logr25))*table[need].d25_maj
          endif

; position angle

          if tag_exist(leda,'pa') then begin
             need = where((table.optical_posangle eq -999.0) and (strcompress(leda.pa,/remove) ne ''),nneed)
             if (nneed ne 0L) then begin
                table[need].optical_posangle = leda[need].pa
                table[need].optical_posangle_origin = 'LEDA'
                table[need].optical_posangle_ref = 'LEDA'
             endif
          endif

; do not use LEDA's inclination angles because it is computed with
; some funny, type-dependent relation
          
;         if tag_exist(leda,'incl') then begin
;            need = where((table.inclination eq -999.0) and (strcompress(leda.incl,/remove) ne ''),nneed)
;            if (nneed ne 0L) then table[need].inclination = leda[need].incl
;         endif

       endif    
       
; other optical diameters: ESO, UGC, MCG, NED "Basic Data"

; ESO    

       good = where((table.d25_maj lt -900.0) and (diamdata.eso_major_axis gt -900.0) and $
         (diamdata.eso_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].d25_maj = diamdata[good].eso_major_axis/60.0
          table[good].d25_min = diamdata[good].eso_minor_axis/60.0
          table[good].d25_origin = 'ESO'
          table[good].d25_ref = diamdata[good].eso_reference
;         struct_print, struct_trimtags(diamdata[good],select=['GALAXY','ESO*'])
       endif

       good = where((table.optical_posangle lt -900.0) and $
         (diamdata.eso_posangle gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].optical_posangle = diamdata[good].eso_posangle
          table[good].optical_posangle_origin = 'ESO'
          table[good].optical_posangle_ref = diamdata[good].eso_reference
       endif

; UGC    
       
       good = where((table.d25_maj lt -900.0) and (diamdata.ugc_major_axis gt -900.0) and $
         (diamdata.ugc_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].d25_maj = diamdata[good].ugc_major_axis/60.0
          table[good].d25_min = diamdata[good].ugc_minor_axis/60.0
          table[good].d25_origin = 'UGC'
          table[good].d25_ref = diamdata[good].ugc_reference
;         struct_print, struct_trimtags(diamdata[good],select=['GALAXY','UGC*'])
       endif

       good = where((table.optical_posangle lt -900.0) and (diamdata.ugc_posangle gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].optical_posangle = diamdata[good].ugc_posangle
          table[good].optical_posangle_origin = 'UGC'
          table[good].optical_posangle_ref = diamdata[good].ugc_reference
       endif

; MCG    
       
       good = where((table.d25_maj lt -900.0) and (diamdata.mcg_major_axis gt -900.0) and $
         (diamdata.mcg_minor_axis gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].d25_maj = diamdata[good].mcg_major_axis/60.0
          table[good].d25_min = diamdata[good].mcg_minor_axis/60.0
          table[good].d25_origin = 'MCG'
          table[good].d25_ref = diamdata[good].mcg_reference
;         struct_print, struct_trimtags(diamdata[good],select=['GALAXY','MCG*'])
       endif

       good = where((table.optical_posangle lt -900.0) and (diamdata.mcg_posangle gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].optical_posangle = diamdata[good].mcg_posangle
          table[good].optical_posangle_origin = 'MCG'
          table[good].optical_posangle_ref = diamdata[good].mcg_reference
       endif

; NED "Basic Data"    
       
       good = where((table.d25_maj lt -900.0) and (basicdata.dmaj gt -900.0) and $
         (basicdata.dmin gt -900.0),ngood)
       if (ngood ne 0L) then begin
          table[good].d25_maj = basicdata[good].dmaj
          table[good].d25_min = basicdata[good].dmin
          table[good].d25_origin = 'NED'
          table[good].d25_ref = 'NED'
;         struct_print, struct_trimtags(basicdata[good])
       endif

; compute the inclination angles

       good = where((table.twomass_maj gt -900.0) and (table.twomass_min gt -900.0),ngood)
       if (ngood ne 0L) then begin
          ratio = table[good].twomass_min/table[good].twomass_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)
          table[good].twomass_inclination = asin(quantity<1)*!radeg
       endif
       
       good = where((table.d25_maj gt -900.0) and (table.d25_min gt -900.0),ngood)
       if (ngood ne 0L) then begin
          ratio = table[good].d25_min/table[good].d25_maj ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2)/0.96)
          table[good].optical_inclination = asin(quantity<1)*!radeg
       endif
       
; now merge the optical and infrared position and inclination angles;
; I found that in general the optical position angles are good, even
; though they are based on "by-eye" measurements

       table.inclination     = table.optical_inclination
       table.posangle        = table.optical_posangle
       table.posangle_origin = table.optical_posangle_origin
       table.posangle_ref    = table.optical_posangle_ref
;      table.posangle        = table.twomass_posangle

       more = where((table.inclination lt -900.0) and (table.twomass_inclination gt -900.0),nmore)
       if (nmore ne 0L) then $
         table[more].inclination = table[more].twomass_inclination

       more = where((table.posangle lt -900.0) and (table.twomass_posangle gt -900.0),nmore)
       if (nmore ne 0L) then begin
          table[more].posangle = table[more].twomass_posangle
          table[more].posangle_origin = '2MASS/'+strtrim(table[more].twomass_jhk_ref,2)
          table[more].posangle_ref = table[more].twomass_diameter_ref
       endif
;      more = where((table.posangle lt -900.0) and (table.optical_posangle gt -900.0),nmore)
;      if (nmore ne 0L) then $
;        table[more].posangle = table[more].optical_inclination

; print messages    
       
       need = where(table.d25_maj lt -900.0,nneed)
       if (nneed ne 0L) then begin
          splog, 'The following objects need optical diameter measurements.'
          niceprint, table[need].galaxy
       endif
       
       need = where(table.posangle lt -900.0,nneed)
       if (nneed ne 0L) then begin
          splog, 'The following objects need position angles.'
          niceprint, table[need].galaxy
       endif
       
       need = where(table.inclination lt -900.0,nneed)
       if (nneed ne 0L) then begin
          splog, 'The following objects need inclination angles.'
          niceprint, table[need].galaxy
       endif

    endif
       
;   splog, 'The following objects need diameters:'
;   need = where(table.d25_maj lt -900.0,nneed)
;   if (nneed ne 0L) then niceprint, table[need].galaxy, table[need].ned_galaxy
    
; compute the FIR flux and luminosity using the Helou et al. (1988)
; formula, which only needs 60 and 100 micron IRAS fluxes.  compute
; L(IR) using two methods: direct integration and the bolometric
; correction technique by Dale et al. (2001); include a 30%
; uncertainty in this bolometric correction; also estimate the
; infrared star-formation rate according to Kennicutt (1998); to
; convert to cgs units: 1 erg/s/cm2 = 1D3 W/m2

    splog, 'Computing IR/FIR fluxes and luminosities.'
    
    good = where((table.iras[2] gt 0.0) and (table.iras[3] gt 0.0),ngood)
    if (ngood ne 0L) then begin

       table[good].fir_flux_sm96[0] = 1D3*1.26D-14*(2.58*table[good].iras[2]+table[good].iras[3]) ; [erg/s/cm2]
       table[good].fir_flux_sm96[1] = 1D3*1.26D-14*sqrt( (2.58*table[good].iras_err[2])^2 + $  ; [erg/s/cm2]
         table[good].iras_err[3]^2 )

; Dale et al. (2001): bolometric (3-1100 micron) IR luminosity;
; include a 30% uncertainty in the correction

       dale = [0.2738D,-0.0282,0.7281,0.6208,0.9118] ; coefficients from Dale et al. 2001

       xratio = alog10(table[good].iras[2]/table[good].iras[3])
       bolcor = 10^(dale[0] + dale[1]*xratio + dale[2]*xratio^2 + dale[3]*xratio^3 + dale[4]*xratio^4)
       bolcor_err = 0.30*bolcor
       
       table[good].ir_flux_dale01[0] = bolcor * table[good].fir_flux_sm96[0]                     ; [erg/s/cm2]
       table[good].ir_flux_dale01[1] = sqrt( (bolcor*table[good].fir_flux_sm96[1])^2 + $ ; [erg/s/cm2]
         (bolcor_err*table[good].fir_flux_sm96[0])^2 )

    endif

; direct integration; assume a fixed 30% uncertainty

    splog, 'Integrating the IR flux numerically.'
    
    good = where((table.fir_flux_sm96[0] gt -900.0) and (table.iras[2] gt 0.0) and $
      (table.iras[3] gt 0.0) and ((table.iras[0] gt 0.0) or (table.iras_12_predict gt -900.0)) and $
      ((table.iras[1] gt 0.0) or (table.iras_25_predict gt -900.0)),ngood)

    if (ngood ne 0L) then begin

       f60 = table[good].iras[2]
       f60_err = table[good].iras_err[2]

       f100 = table[good].iras[3]
       f100_err = table[good].iras_err[3]

       f12 = f60*0.0
       f12_err = f60*0.0

       predict = where((table[good].iras_12_predict gt -900) and $
         (table[good].iras[0] lt 0.0),npredict,comp=measured,ncomp=nmeasured)
       if (npredict ne 0L) then begin
          f12[predict] = table[good[predict]].iras_12_predict
          f12_err[predict] = table[good[predict]].iras_12_predict_err
       endif
       if (nmeasured ne 0L) then begin
          f12[measured] = table[good[measured]].iras[0]
          f12_err[measured] = table[good[measured]].iras_err[0]
       endif
       
       f25 = f60*0.0
       f25_err = f60*0.0

       predict = where((table[good].iras_25_predict gt -900) and $
         (table[good].iras[1] lt 0.0),npredict,comp=measured,ncomp=nmeasured)
       if (npredict ne 0L) then begin
          f25[predict] = table[good[predict]].iras_25_predict
          f25_err[predict] = table[good[predict]].iras_25_predict_err
       endif
       if (nmeasured ne 0L) then begin
          f25[measured] = table[good[measured]].iras[1]
          f25_err[measured] = table[good[measured]].iras_err[1]
       endif

       iras_flux = transpose([ [f12], [f25], [f60], [f100] ])
       iras_ferr = transpose([ [f12_err], [f25_err], [f60_err], [f100_err] ])

       result = im_total_ir(iras_flux,iras_ferr=iras_ferr,alpha=alpha)
;      result = im_total_ir(iras_flux[*,0],iras_ferr=iras_ferr[*,0],alpha=alpha)
       
       table[good].ir_flux[0]   = result.ir_flux
       table[good].ir_flux[1]   = result.ir_ferr
       table[good].ir_dust_temp = result.dust_temp
       
    endif
       
; bolometric IR flux from Dale & Helou (2002)
    
    good = where((table.iras[1] gt 0.0) and (table.iras[2] gt 0.0) and $
      (table.iras[3] gt 0.0) and (table.fir_flux_sm96[0] gt -900),ngood)
    if (ngood ne 0L) then begin

       nu = light*1D9/[25.0,60.0,100.0]
       dale = [2.403D,-0.2454D,1.6381D] ; coefficients from Dale & Helou 2002

       table[good].ir_flux_dale02[0] = (dale[0]*nu[0]*table[good].iras[1] + $
         dale[1]*nu[1]*table[good].iras[2] + dale[2]*nu[2]*table[good].iras[3])*1D-23
       table[good].ir_flux_dale02[1] = sqrt(((dale[0]*nu[0]*table[good].iras_err[1])^2 + $
         (dale[1]*nu[1]*table[good].iras_err[2])^2 + (dale[2]*nu[2]*table[good].iras_err[3])^2))*1D-23
       
    endif

; bolometric IR flux from Sanders & Mirabel (1996)
    
    good = where((table.iras[0] gt 0.0) and (table.iras[1] gt 0.0) and $
      (table.iras[2] gt 0.0) and (table.iras[3] gt 0.0) and $
      (table.fir_flux_sm96[0] gt -900),ngood)
    if (ngood ne 0L) then begin

       table[good].ir_flux_sm96[0] = 1D3*1.8D-14*(13.48*table[good].iras[0]+$
         5.16*table[good].iras[1]+2.58*table[good].iras[2]+table[good].iras[3]) ; [erg/s/cm2]
       table[good].ir_flux_sm96[1] = 1D3*1.8D-14*sqrt(((13.48*table[good].iras_err[0])^2.0+$
         (5.16*table[good].iras_err[1])^2.0+(2.58*table[good].iras_err[2])^2.0+table[good].iras_err[3]^2)) ; [erg/s/cm2]
       
    endif

; compute the IR and FIR luminosity and error

    good = where((table.distance gt -900.0) and (table.fir_flux_sm96[0] gt -900.0),ngood)
    if (ngood ne 0L) then begin
       fir_lum = table[good].fir_flux_sm96[0]*4.0*!dpi*(table[good].distance*mpc2cm)^2.0                         ; [erg/s]
       fir_lum_err = 4.0*!dpi * sqrt( (table[good].fir_flux_sm96[1]*(table[good].distance*mpc2cm)^2)^2 + $
         (table[good].fir_flux_sm96[0]*(2*table[good].distance*mpc2cm*table[good].distance_err*mpc2cm))^2 ) ; [erg/s]

       table[good].fir_lum_sm96[1] = fir_lum_err/fir_lum/alog(10.0) ; log L_sun
       table[good].fir_lum_sm96[0] = alog10(fir_lum/lsun)               ; log L_sun
    endif

    good = where((table.distance gt -900.0) and (table.ir_flux[0] gt -900.0),ngood)
    if (ngood ne 0L) then begin
       ir_lum = table[good].ir_flux[0]*4.0*!dpi*(table[good].distance*mpc2cm)^2.0                           ; [erg/s]
       ir_lum_err = 4.0*!dpi * sqrt( (table[good].ir_flux[1]*(table[good].distance*mpc2cm)^2)^2 + $
         (table[good].ir_flux[0]*(2*table[good].distance*mpc2cm*table[good].distance_err*mpc2cm))^2 )  ; [erg/s]

       table[good].ir_lum[1] = ir_lum_err/ir_lum/alog(10.0) ; log L_sun
       table[good].ir_lum[0] = alog10(ir_lum/lsun)              ; log L_sun
    endif

    good = where((table.distance gt -900.0) and (table.ir_flux_dale02[0] gt -900.0),ngood)
    if (ngood ne 0L) then begin
       ir_lum_dale02 = table[good].ir_flux_dale02[0]*4.0*!dpi*(table[good].distance*mpc2cm)^2.0                           ; [erg/s]
       ir_lum_dale02_err = 4.0*!dpi * sqrt( (table[good].ir_flux_dale02[1]*(table[good].distance*mpc2cm)^2)^2 + $
         (table[good].ir_flux_dale02[0]*(2*table[good].distance*mpc2cm*table[good].distance_err*mpc2cm))^2 )  ; [erg/s]

       table[good].ir_lum_dale02[1] = ir_lum_dale02_err/ir_lum_dale02/alog(10.0) ; log L_sun
       table[good].ir_lum_dale02[0] = alog10(ir_lum_dale02/lsun)              ; log L_sun
    endif

    good = where((table.distance gt -900.0) and (table.ir_flux_dale01[0] gt -900.0),ngood)
    if (ngood ne 0L) then begin
       ir_lum_dale01 = table[good].ir_flux_dale01[0]*4.0*!dpi*(table[good].distance*mpc2cm)^2.0                           ; [erg/s]
       ir_lum_dale01_err = 4.0*!dpi * sqrt( (table[good].ir_flux_dale01[1]*(table[good].distance*mpc2cm)^2)^2 + $
         (table[good].ir_flux_dale01[0]*(2*table[good].distance*mpc2cm*table[good].distance_err*mpc2cm))^2 )  ; [erg/s]

       table[good].ir_lum_dale01[1] = ir_lum_dale01_err/ir_lum_dale01/alog(10.0) ; log L_sun
       table[good].ir_lum_dale01[0] = alog10(ir_lum_dale01/lsun)              ; log L_sun
    endif

    good = where((table.distance gt -900.0) and (table.ir_flux_sm96[0] gt -900.0),ngood)
    if (ngood ne 0L) then begin
       ir_lum_sm96 = table[good].ir_flux_sm96[0]*4.0*!dpi*(table[good].distance*mpc2cm)^2.0                           ; [erg/s]
       ir_lum_sm96_err = 4.0*!dpi * sqrt( (table[good].ir_flux_sm96[1]*(table[good].distance*mpc2cm)^2)^2 + $
         (table[good].ir_flux_sm96[0]*(2*table[good].distance*mpc2cm*table[good].distance_err*mpc2cm))^2 )  ; [erg/s]

       table[good].ir_lum_sm96[1] = ir_lum_sm96_err/ir_lum_sm96/alog(10.0) ; log L_sun
       table[good].ir_lum_sm96[0] = alog10(ir_lum_sm96/lsun)              ; log L_sun
    endif

; sort by RA

    if (not keyword_set(norasort)) then begin
       splog, 'Sorting the output by right ascension.'
       srtra = sort(im_hms2dec(table.ra))
       table = table[srtra]
    endif

; write out the full atlas table as a binary fits file if requested 

    if keyword_set(write) then begin
       splog, 'Writing '+outpath+outname+'.'
       mwrfits, table, outpath+outname, /create
       spawn, ['gzip -f '+outpath+outname], /sh
    endif

return
end

; ###########################################################################
; OLD CODE THAT MAY BE USEFUL LATER
; ###########################################################################

; based on the literature, determine which objects have nuclear AGN;
; RELEGATED: use these samples on a case-by-case basis (jm08jun06nyu)

;;; read the Ho et al (1997) database    
;;
;;    ho97 = read_97ho()
;;    raho = 15.0D*im_hms2dec(ho97.ra)
;;    decho = im_hms2dec(ho97.dec)
;;    
;;    ntot = im_djs_angle_match(ra,dec,raho,decho,dtheta=1.0/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;;    match = where(mindx ne -1L,nmatch)
;;    if nmatch ne 0L then begin
;;       table[match].lit_class = strcompress(ho97[mindx[match]].class,/remove)
;;       table[match].lit_class_ref = 'Ho et al. (1997)'
;;    endif
;;
;;; read the Veilleux et al (1995, 1999) database    
;;
;;    veilleux = read_99veilleux()
;;    raveilleux = 15.0D*im_hms2dec(veilleux.ra)
;;    decveilleux = im_hms2dec(veilleux.dec)
;;    
;;    ntot = im_djs_angle_match(ra,dec,raveilleux,decveilleux,dtheta=1.0/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;;    match = where(mindx ne -1L,nmatch)
;;    if nmatch ne 0L then begin
;;       for i = 0L, nmatch-1L do begin
;;          if table[match[i]].lit_class eq '' then begin
;;             table[match[i]].lit_class = strtrim(veilleux[mindx[match[i]]].class,2)
;;             table[match[i]].lit_class_ref = 'Veilleux et al. (1995,1999)'
;;          endif else begin
;;             table[match[i]].lit_class = table[match[i]].lit_class+', '+strtrim(veilleux[mindx[match[i]]].class,2)
;;             table[match[i]].lit_class_ref = table[match[i]].lit_class_ref+', Veilleux et al. (1995,1999)'
;;          endelse
;;       endfor
;;    endif
;;
;;;   struct_print, struct_trimtags(table,select=['galaxy','lit_class','lit_class_ref'])
;;
;;    nolit = where(strcompress(table.lit_class,/remove) eq '',nnolit)
;;;   splog, 'There are '+string(nnolit,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
;;;     ' galaxies with no literature classification.'
    

;;; ##########
;;; IRAS 12
;;; ##########
;;
;;       good = where((table.iras_12 lt -900) and (photodata.sanders_iras_12 gt -900.0) and $
;;         (photodata.sanders_iras_12_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_12     = photodata[good].sanders_iras_12
;;          table[good].iras_12_err = photodata[good].sanders_iras_12_err
;;          table[good].iras_12_ref = 'Sanders et al. 2003'
;;       endif
;;       
;;       good = where((table.iras_12 lt -900) and (photodata.rice_iras_12 gt -900.0) and $
;;         (photodata.rice_iras_12_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_12     = photodata[good].rice_iras_12
;;          table[good].iras_12_err = photodata[good].rice_iras_12_err
;;          table[good].iras_12_ref = 'Rice et al. 1988'
;;       endif
;;       
;;       good = where((table.iras_12 lt -900) and (photodata.soifer_iras_12 gt -900.0) and $
;;         (photodata.soifer_iras_12_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_12     = photodata[good].soifer_iras_12
;;          table[good].iras_12_err = photodata[good].soifer_iras_12_err
;;          table[good].iras_12_ref = 'Soifer et al. 1989'
;;       endif
;;
;;       good = where((table.iras_12 lt -900) and (photodata.moshir_iras_12 gt -900.0) and $
;;         (photodata.moshir_iras_12_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_12     = photodata[good].moshir_iras_12
;;          table[good].iras_12_err = photodata[good].moshir_iras_12_err
;;          table[good].iras_12_ref = 'Moshir et al. 1990'
;;       endif
;;       
;;; ##########
;;; IRAS 25
;;; ##########
;;
;;       good = where((table.iras_25 lt -900) and (photodata.sanders_iras_25 gt -900.0) and $
;;         (photodata.sanders_iras_25_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_25     = photodata[good].sanders_iras_25
;;          table[good].iras_25_err = photodata[good].sanders_iras_25_err
;;          table[good].iras_25_ref = 'Sanders et al. 2003'
;;       endif
;;       
;;       good = where((table.iras_25 lt -900) and (photodata.rice_iras_25 gt -900.0) and $
;;         (photodata.rice_iras_25_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_25     = photodata[good].rice_iras_25
;;          table[good].iras_25_err = photodata[good].rice_iras_25_err
;;          table[good].iras_25_ref = 'Rice et al. 1988'
;;       endif
;;       
;;       good = where((table.iras_25 lt -900) and (photodata.soifer_iras_25 gt -900.0) and $
;;         (photodata.soifer_iras_25_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_25     = photodata[good].soifer_iras_25
;;          table[good].iras_25_err = photodata[good].soifer_iras_25_err
;;          table[good].iras_25_ref = 'Soifer et al. 1989'
;;       endif
;;
;;       good = where((table.iras_25 lt -900) and (photodata.moshir_iras_25 gt -900.0) and $
;;         (photodata.moshir_iras_25_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_25     = photodata[good].moshir_iras_25
;;          table[good].iras_25_err = photodata[good].moshir_iras_25_err
;;          table[good].iras_25_ref = 'Moshir et al. 1990'
;;       endif
;;       
;;; ##########
;;; IRAS 60
;;; ##########
;;
;;       good = where((table.iras_60 lt -900) and (photodata.sanders_iras_60 gt -900.0) and $
;;         (photodata.sanders_iras_60_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_60     = photodata[good].sanders_iras_60
;;          table[good].iras_60_err = photodata[good].sanders_iras_60_err
;;          table[good].iras_60_ref = 'Sanders et al. 2003'
;;       endif
;;       
;;       good = where((table.iras_60 lt -900) and (photodata.rice_iras_60 gt -900.0) and $
;;         (photodata.rice_iras_60_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_60     = photodata[good].rice_iras_60
;;          table[good].iras_60_err = photodata[good].rice_iras_60_err
;;          table[good].iras_60_ref = 'Rice et al. 1988'
;;       endif
;;       
;;       good = where((table.iras_60 lt -900) and (photodata.soifer_iras_60 gt -900.0) and $
;;         (photodata.soifer_iras_60_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_60     = photodata[good].soifer_iras_60
;;          table[good].iras_60_err = photodata[good].soifer_iras_60_err
;;          table[good].iras_60_ref = 'Soifer et al. 1989'
;;       endif
;;
;;       good = where((table.iras_60 lt -900) and (photodata.moshir_iras_60 gt -900.0) and $
;;         (photodata.moshir_iras_60_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_60     = photodata[good].moshir_iras_60
;;          table[good].iras_60_err = photodata[good].moshir_iras_60_err
;;          table[good].iras_60_ref = 'Moshir et al. 1990'
;;       endif
;;       
;;; ##########
;;; IRAS 100
;;; ##########
;;
;;       good = where((table.iras_100 lt -900) and (photodata.sanders_iras_100 gt -900.0) and $
;;         (photodata.sanders_iras_100_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_100     = photodata[good].sanders_iras_100
;;          table[good].iras_100_err = photodata[good].sanders_iras_100_err
;;          table[good].iras_100_ref = 'Sanders et al. 2003'
;;       endif
;;       
;;       good = where((table.iras_100 lt -900) and (photodata.rice_iras_100 gt -900.0) and $
;;         (photodata.rice_iras_100_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_100     = photodata[good].rice_iras_100
;;          table[good].iras_100_err = photodata[good].rice_iras_100_err
;;          table[good].iras_100_ref = 'Rice et al. 1988'
;;       endif
;;       
;;       good = where((table.iras_100 lt -900) and (photodata.soifer_iras_100 gt -900.0) and $
;;         (photodata.soifer_iras_100_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_100     = photodata[good].soifer_iras_100
;;          table[good].iras_100_err = photodata[good].soifer_iras_100_err
;;          table[good].iras_100_ref = 'Soifer et al. 1989'
;;       endif
;;
;;       good = where((table.iras_100 lt -900) and (photodata.moshir_iras_100 gt -900.0) and $
;;         (photodata.moshir_iras_100_err gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;          table[good].iras_100     = photodata[good].moshir_iras_100
;;          table[good].iras_100_err = photodata[good].moshir_iras_100_err
;;          table[good].iras_100_ref = 'Moshir et al. 1990'
;;       endif

; the following exception is obsolete with the introduction of the
; Sanders et al. 2003 IRAS fluxes (jm07nov29nyu)
    
;; make an exception for NGC5194 (M51a); for this object, the Soifer
;; fluxes are closer to the truth (Kennicutt 2005, private
;; communication) and the Moshir fluxes are off by an
;; order-of-magnitude (especially at 12 and 25 microns); however,
;; Soifer does not give a flux error, so assume 25%
;
;    indx = where(strmatch(strcompress(table.galaxy,/remove),'*ngc5194*',/fold) eq 1B,nindx)
;    if (nindx ne 0L) then begin
;       table[indx].iras_12      = photodata[indx].soifer_iras_12
;       table[indx].iras_12_err  = photodata[indx].soifer_iras_12*0.25
;       table[indx].iras_12_ref  = 'Soifer et al. 1989'
;       table[indx].iras_25      = photodata[indx].soifer_iras_25
;       table[indx].iras_25_err  = photodata[indx].soifer_iras_25*0.25
;       table[indx].iras_25_ref  = 'Soifer et al. 1989'
;       table[indx].iras_60      = photodata[indx].soifer_iras_60
;       table[indx].iras_60_err  = photodata[indx].soifer_iras_60*0.25
;       table[indx].iras_60_ref  = 'Soifer et al. 1989'
;       table[indx].iras_100     = photodata[indx].soifer_iras_100
;       table[indx].iras_100_err = photodata[indx].soifer_iras_100*0.25
;       table[indx].iras_100_ref = 'Soifer et al. 1989'
;    endif
    
;;; literature gas masses; REMOVE (jm08june06nyu)
;;
;;    splog, 'Compiling gas masses from the literature.'
;;
;;    davoust04 = read_04davoust()
;;    radav = 15.0D*im_hms2dec(davoust04._raj2000)
;;    decdav = im_hms2dec(davoust04._dej2000)
;;    
;;    ntot = im_djs_angle_match(ra,dec,radav,decdav,dtheta=10.0/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;;    match = where(mindx ne -1L,nmatch)
;;;   niceprint, table[match].galaxy, table[match].ned_galaxy, $
;;;     'MRK'+string(davoust04[mindx[match]].mrk,format='(I0)')
;;
;;; correct their masses to our distances
;;    
;;    if (nmatch ne 0L) then begin
;;       good = where(table[match].distance gt -900,ngood)
;;       if (ngood ne 0L) then begin
;;          table[match[good]].mass_HI_lit_err = davoust04[mindx[match[good]]].e_logMHI
;;          table[match[good]].mass_HI_lit = davoust04[mindx[match[good]]].logMHI + $
;;            2*alog10(table[match[good]].distance/davoust04[mindx[match[good]]].dist)
;;          table[match[good]].lit_class_ref = 'Davoust et al. 2004'
;;       endif
;;    endif
    
;;; additional measurements from LEDA
;;
;;    if (nleda ne 0L) then begin
;;
;;       splog, 'Parsing additional measurements from LEDA'
;;    
;;; hydrogen line widths    
;;       
;;       if tag_exist(leda,'w20') then begin
;;          need = where(strcompress(leda.w20,/remove) ne '',nneed)
;;          if (nneed ne 0L) then table[need].w20 = leda[need].w20
;;       endif
;;
;;       if tag_exist(leda,'e_w20') then begin
;;          need = where(strcompress(leda.e_w20,/remove) ne '',nneed)
;;          if (nneed ne 0L) then table[need].w20_err = leda[need].e_w20
;;       endif
;;
;;;;    if tag_exist(leda,'w50') then begin
;;;;       need = where(strcompress(leda.w50,/remove) ne '',nneed)
;;;;       if (nneed ne 0L) then table[need].w50 = leda[need].w50
;;;;    endif
;;
;;; compute the hydrogen mass, corrected for self-absorption
;;
;;       if tag_exist(leda,'m21c') and tag_exist(leda,'e_m21') then begin
;;
;;          need = where((strcompress(leda.m21c,/remove) ne '') and (strcompress(leda.e_m21,/remove) ne '') and $
;;            (table.distance ne -999.0) and (table.distance_err ne -999.0),nneed)
;;          if (nneed ne 0L) then begin
;;
;;             f21 = 10.0^(-0.4*(leda[need].m21c-17.40)) ; [Jy km/s]
;;             sf21 = 0.4 * alog(10.0) * f21 * leda[need].e_m21
;;
;;             d = table[need].distance
;;             d_err = table[need].distance_err
;;             MHI = 2.36E5 * d^2.0 * f21 ; [M_sun]
;;             MHI_err = 2.36E5 * sqrt( (2*d*f21*d_err)^2.0 + ((d^2.0*sf21)^2.0) )
;;             
;;             table[need].mass_HI_leda_err = MHI_err/MHI/alog(10.0) ; log M_HI/M_sun
;;             table[need].mass_HI_leda = alog10(MHI)
;;
;;          endif 
;;       endif
;;
;;    endif

;; compute the gas fraction using relations in McGaugh & de Block
;; (1997) 
;
;    good = where(table.rc3_t gt -900.0,ngood)
;    if (ngood ne 0L) then begin
;       type = table[good].rc3_t
;       gas_ratio = ((3.7-0.8*type+0.043*type^2)>0.0)<4.0 ; gas_ratio --> 0 as type --> 10
;       table[good].eta_gas = 1.4*(1+gas_ratio)
;    endif
;
;    good = where((table.eta_gas gt -900.0) and (table.mass_HI gt -900.0),ngood)
;    if (ngood ne 0L) then begin
;       table[good].mass_gas     = table[good].mass_HI + alog10(table[good].eta_gas)
;       table[good].mass_gas_err = table[good].mass_HI_err
;    endif
    
;;; FUV photometry, radio data, errors, and associated quantities 
;;
;;    splog, 'Compiling FUV and radio photometry.'
;;    
;;; cross-match with the Bell (2003) catalog to derive FUV and 1.4 GHz
;;; data for our sample
;;
;;    bell = read_03bell_radio_fir()
;;    
;;    bellra = 15.0D*im_hms2dec(bell.ra)
;;    belldec = im_hms2dec(bell.dec)
;;
;;    splog, 'Bell (2002) Radio-FIR search radius = '+string(1.0,format='(G0.0)')+'".'
;;    ntot = im_djs_angle_match(ra,dec,bellra,belldec,dtheta=1.0/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;;    
;;    match = where(mindx ne -1L,nmatch)
;;    if (nmatch ne 0L) then begin
;;
;;       splog, 'Found '+string(nmatch,format='(I0)')+' objects in common.'
;;;      niceprint, table[match].ned_galaxy, bell[mindx[match]].ned_galaxy
;;
;;       dist = bell[mindx[match]].bell_dis*mpc2cm ; [cm]
;;       area = 4.0*!dpi*dist*dist                 ; [cm2]
;;
;;; FUV fluxes and luminosities near ~1550 Angstrom; the 1D7 converts
;;; W to erg/s
;;       
;;       fuvgood = where(bell[mindx[match]].bell_fuv gt -900.0,nfuvgood)
;;       if (nfuvgood ne 0L) then begin
;;
;;          table[match[fuvgood]].fuv_ref = 'Bell (2003)'
;;          
;;          fuvlum = 1550.0*10^bell[mindx[match[fuvgood]]].bell_fuv               ; [W]
;;          fuvlum_err = bell[mindx[match[fuvgood]]].bell_e_fuv*fuvlum*alog(10.0) ; [W]
;;
;;          table[match[fuvgood]].fuv_flux = fuvlum*1D7/area[fuvgood]         ; [erg/s/cm2]
;;          table[match[fuvgood]].fuv_flux_err = fuvlum_err*1D7/area[fuvgood] ; [erg/s/cm2]
;;
;;          good = where(table[match[fuvgood]].distance gt -900.0,ngood)
;;          if (ngood ne 0L) then begin
;;
;;             fuv_lum = table[match[fuvgood[good]]].fuv_flux * $            
;;               4.0*!dpi*(table[match[fuvgood[good]]].distance*mpc2cm)^2.0                 ; [erg/s]
;;             
;;             fuv_lum_err = 4.0*!dpi * sqrt( (table[match[fuvgood[good]]].fuv_flux_err * $ ; [erg/s]
;;               (table[match[fuvgood[good]]].distance*mpc2cm)^2)^2 + $
;;               (table[match[fuvgood[good]]].fuv_flux*(2*table[match[fuvgood[good]]].distance * $
;;               mpc2cm*table[match[fuvgood[good]]].distance_err*mpc2cm))^2 ) 
;;             
;;; take the log of the quantities
;;             
;;             table[match[fuvgood[good]]].fuv_lum_err = fuv_lum_err/fuv_lum/alog(10.0) ; log [L_sun]
;;             table[match[fuvgood[good]]].fuv_lum = alog10(fuv_lum/lsun)               ; log [L_sun]
;;
;;          endif 
;;          
;;       endif
;;
;;; radio 1.4GHz fluxes and luminosities
;;       
;;       radiogood = where(bell[mindx[match]].bell_1_4GHZ gt -900.0,nradiogood)
;;       if (nradiogood ne 0L) then begin
;;          
;;          radiolum = 1.4D9*10^bell[mindx[match[radiogood]]].bell_1_4GHZ                  ; [W]
;;          radiolum_err = bell[mindx[match[radiogood]]].bell_e_1_4GHZ*radiolum*alog(10.0) ; [W]
;;
;;          table[match[radiogood]].radio_1_4GHZ_flux = radiolum*1D7/area[radiogood]         ; [erg/s/cm2]
;;          table[match[radiogood]].radio_1_4GHZ_flux_err = radiolum_err*1D7/area[radiogood] ; [erg/s/cm2]
;;
;;          good = where(table[match[radiogood]].distance gt -900.0,ngood)
;;          if (ngood ne 0L) then begin
;;
;;             radio_lum = table[match[radiogood[good]]].radio_1_4GHZ_flux * $            
;;               4.0*!dpi*(table[match[radiogood[good]]].distance*mpc2cm)^2.0                            ; [erg/s]
;;             
;;             radio_lum_err = 4.0*!dpi * sqrt( (table[match[radiogood[good]]].radio_1_4GHZ_flux_err * $ ; [erg/s]
;;               (table[match[radiogood[good]]].distance*mpc2cm)^2)^2 + $
;;               (table[match[radiogood[good]]].radio_1_4GHZ_flux*(2*table[match[radiogood[good]]].distance * $
;;               mpc2cm*table[match[radiogood[good]]].distance_err*mpc2cm))^2 ) 
;;             
;;; take the log of the quantities
;;             
;;             table[match[radiogood[good]]].radio_1_4GHZ_lum_err = radio_lum_err/radio_lum/alog(10.0) ; log [L_sun]
;;             table[match[radiogood[good]]].radio_1_4GHZ_lum = alog10(radio_lum/lsun)                 ; log [L_sun]
;;
;;          endif 
;;       endif
;;       
;;    endif 

;;; supplement the Bell (2003) FUV fluxes with data from Rifatto et
;;; al. (1995)
;;
;;    good = where((table.fuv_flux lt -900) and (photodata.uv_1650 gt -900),ngood)
;;    if (ngood ne 0L) then begin
;;
;;; reject galaxies that are too large; assume that galaxies without
;;; measured diameters are small enough; the area rejection criterion is
;;; from Bell (2003)
;;       
;;       Dgood = where((table[good].d25_maj gt -900) and (table[good].d25_min gt -900),nDgood)
;;       if (nDgood ne 0L) then begin
;;          keep = where(alog10(table[good[Dgood]].d25_maj*table[good[Dgood]].d25_min+1) lt 1.5,nkeep)
;;          if (nkeep ne 0L) then begin
;;             table[good[Dgood[keep]]].fuv_flux     = 1650.0*photodata[good[Dgood[keep]]].uv_1650     ; [erg/s/cm2]
;;             table[good[Dgood[keep]]].fuv_flux_err = 1650.0*photodata[good[Dgood[keep]]].uv_1650_err ; [erg/s/cm2]
;;             table[good[Dgood[keep]]].fuv_ref = 'Rifatto et al. (1995)'
;;          endif
;;          
;;       endif
;;       
;;       small = where((table[good].d25_maj lt -900) or (table[good].d25_min lt -900),nsmall)
;;       if (nsmall ne 0L) then begin
;;          table[good[small]].fuv_flux     = 1650.0*photodata[good[small]].uv_1650     ; [erg/s/cm2]
;;          table[good[small]].fuv_flux_err = 1650.0*photodata[good[small]].uv_1650_err ; [erg/s/cm2]
;;          table[good[small]].fuv_ref = 'Rifatto et al. (1995)'
;;       endif
;;
;;; compute luminosities
;;
;;       gdist = where((table[good].distance gt -900.0) and (table[good].fuv_flux gt -900),ngdist)
;;       if (ngdist ne 0L) then begin
;;
;;          fuv_lum = table[good[gdist]].fuv_flux * $            
;;            4.0*!dpi*(table[good[gdist]].distance*mpc2cm)^2.0                 ; [erg/s]
;;          
;;          fuv_lum_err = 4.0*!dpi * sqrt( (table[good[gdist]].fuv_flux_err * $ ; [erg/s/A]
;;            (table[good[gdist]].distance*mpc2cm)^2)^2 + $
;;            (table[good[gdist]].fuv_flux*(2*table[good[gdist]].distance * $
;;            mpc2cm*table[good[gdist]].distance_err*mpc2cm))^2 ) 
;;          
;;; take the log of the quantities
;;          
;;          table[good[gdist]].fuv_lum_err = fuv_lum_err/fuv_lum/alog(10.0) ; log [L_sun]
;;          table[good[gdist]].fuv_lum = alog10(fuv_lum/lsun) ; log [L_sun]
;;
;;       endif 
;;          
;;    endif 

;; ---------------------------------------------------------------------------
;; optical/FIR and optical/IR ratios and errors; since the B-band has
;; been normalized to the B-band luminosity, apply the appropriate
;; conversion 
;; ---------------------------------------------------------------------------
;
;    splog, 'Computing luminosity ratios.'
;    
;    good = where((table.rc3_B_lum gt -900.0) and (table.fir_lum gt -900),ngood)
;    if (ngood ne 0L) then begin
;
;       B_lum = table[good].rc3_B_lum + 0.4*(mbolsun-bmsun)
;
;       table[good].l_fir_l_b = 10^(table[good].fir_lum-B_lum)
;       table[good].l_ir_l_b  = 10^(table[good].ir_lum-B_lum)
;
;       x1 = 10^table[good].fir_lum  & x1err = table[good].fir_lum_err*alog(10.0)*x1
;       y1 = 10^B_lum                & y1err = table[good].rc3_B_lum_err*alog(10.0)*y1
;       table[good].l_fir_l_b_err = im_compute_error(x1,x1err,y1,y1err,/quotient)
;
;       x1 = 10^table[good].ir_lum  & x1err = table[good].ir_lum_err*alog(10.0)*x1
;       y1 = 10^B_lum               & y1err = table[good].rc3_B_lum_err*alog(10.0)*y1
;       table[good].l_ir_l_b_err = im_compute_error(x1,x1err,y1,y1err,/quotient)
;
;    endif

;;; compute the extinction correction to face-on inclination in the
;;; B-band; bootstrap the correction to the U and V-bands using the
;;; O'Donnel Galactic extinction curve
;;
;;    splog, 'Computing the extinction corrections to face-on inclination.'
;;
;;    good = where((table.inclination gt -900.0) and (table.d25_maj gt -900.0) and $
;;      (table.d25_min gt -900.0) and (table.rc3_m_b gt -900.0),ngood)
;;    if (ngood ne 0L) then begin
;;
;;       M_B = table[good].rc3_m_b
;;       ab = table[good].d25_maj/table[good].d25_min ; major-to-minor axis ratio (a/b)
;;       
;;       lolum = where((M_B gt -15.6),nlolum,comp=hilum,ncomp=nhilum)
;;       if (nlolum ne 0L) then begin
;;          table[good[lolum]].A_U_incl = 0.0
;;          table[good[lolum]].A_B_incl = 0.0
;;          table[good[lolum]].A_V_incl = 0.0
;;       endif
;;
;;       if (nhilum ne 0L) then begin
;;
;;          gamma_b = -0.35*(15.6 + M_B[hilum] + 5*alog10(h100/0.8))
;;
;;          table[good[hilum]].A_B_incl = gamma_b*alog10(ab[hilum])
;;          table[good[hilum]].A_U_incl = table[good[hilum]].A_B_incl*k_lambda(U_weff,/odonnell)/k_lambda(B_weff,/odonnell)
;;          table[good[hilum]].A_V_incl = table[good[hilum]].A_B_incl*k_lambda(V_weff,/odonnell)/k_lambda(B_weff,/odonnell)
;;
;;;         niceprint, table[good[hilum]].a_u_incl, table[good[hilum]].a_b_incl, table[good[hilum]].a_v_incl
;;
;;       endif
;;
;;    endif
    
;      table.iras_12     = photodata.iras_12
;      table.iras_12_err = photodata.iras_12_err
;      table.iras_12_ref = photodata.iras_12_ref
;
;      table.iras_25     = photodata.iras_25
;      table.iras_25_err = photodata.iras_25_err
;      table.iras_25_ref = photodata.iras_25_ref
;
;      table.iras_60     = photodata.iras_60
;      table.iras_60_err = photodata.iras_60_err
;      table.iras_60_ref = photodata.iras_60_ref
;
;      table.iras_100     = photodata.iras_100
;      table.iras_100_err = photodata.iras_100_err
;      table.iras_100_ref = photodata.iras_100_ref

;; compute the SFR
;
;;   table[good].sfr_ir = 4.5D-44*ir_lum         ; [M_sun/yr]
;;   table[good].sfr_ir_err = 4.5D-44*ir_lum_err
; 
;;   table[good].sfr_ir_dale01 = 4.5D-44*ir_lum_dale01         ; [M_sun/yr]
;;   table[good].sfr_ir_dale01_err = 4.5D-44*ir_lum_dale01_err

;;; compute UBVJHKs absolute magnitudes and luminosities
;;
;;    splog, 'Computing absolute magnitudes and luminosities.'
;;
;;    band = ['U','B','V','twomass_J','twomass_H','twomass_Ks',$
;;      'twomass_J','twomass_H','twomass_Ks']
;;    
;;    tags = ['rc3_u','rc3_b','rc3_v','twomass_j','twomass_h',$
;;      'twomass_ks','twomass_j20','twomass_h20','twomass_ks20']
;;    tags_err = tags+'_err'
;;
;;    abstags = ['rc3_m_u','rc3_m_b','rc3_m_v','twomass_m_j','twomass_m_h',$
;;      'twomass_m_ks','twomass_m_j20','twomass_m_h20','twomass_m_ks20']
;;    abstags_err = abstags+'_err'
;;    
;;    lumtags = ['rc3_u_lum','rc3_b_lum','rc3_v_lum','twomass_j_lum','twomass_h_lum',$
;;      'twomass_ks_lum','twomass_j20_lum','twomass_h20_lum','twomass_ks20_lum']
;;    lumtags_err = lumtags+'_err'
;;
;;    for iband = 0L, n_elements(tags)-1L do begin
;;
;;       true = tag_exist(table,tags[iband],index=tagsindx)
;;       true = tag_exist(table,tags_err[iband],index=tagsindx_err)
;;       true = tag_exist(table,abstags[iband],index=abstagsindx)
;;       true = tag_exist(table,abstags_err[iband],index=abstagsindx_err)
;;       true = tag_exist(table,lumtags[iband],index=lumtagsindx)
;;       true = tag_exist(table,lumtags_err[iband],index=lumtagsindx_err)
;;
;;       good = where((table.distance gt -900.0) and (table.(tagsindx) gt -900.0),ngood)
;;       if (ngood ne 0L) then begin
;;
;;          mags = im_absolute_magnitudes(band[iband],table[good].(tagsindx),$
;;            table[good].distance,mag_err=table[good].(tagsindx_err),$
;;            distance_err=table[good].distance_err)
;;
;;          table[good].(abstagsindx)     = mags.absmag
;;          table[good].(abstagsindx_err) = mags.absmag_err
;;          
;;          table[good].(lumtagsindx)     = mags.lum
;;          table[good].(lumtagsindx_err) = mags.lum_err
;;          
;;       endif
;;       
;;    endfor

