function init_ancillary_data, ngalaxy=ngalaxy
; jm05may16uofa - excised from WRITE_ANCILLARY_DATA
; jm08apr14nyu - remove some tags that are ~obsolete

    if (n_elements(ngalaxy) eq 0L) then ngalaxy = 1L
    
    table = {$
      galaxy:                          '...', $ ; adopted galaxy name
      ned_galaxy:                      '...', $ ; NED galaxy name
      leda_galaxy:                     '...', $ ; LEDA galaxy name
      ra:                              '...', $ ; RA
      dec:                             '...', $ ; DEC
      position_ref:                    'NED', $ ; position reference
      z:                              -999.0, $ ; redshift and error
      z_err:                          -999.0, $
      cz:                             -999.0, $ ; redshift and error
      cz_err:                         -999.0, $
      xscale:                         -999.0, $ ; [pc/arcsec] at that distance
      ned_morphology:                  '...', $ ; NED morphology and type
      ned_class:                       '...', $ ; NED class
      ned_note:                        '...', $ ; NED special note
      rc3_type:                        '...', $ ; RC3/LEDA morphological type
      rc3_t:                          -999.0, $ ; RC3/LEDA numerical type (continuous)
      lit_type:                        '...', $ ; literature morphological type
      lit_t:                          -999.0, $ ; literature numerical type (discrete)
                              
;     bar:                                0L, $ ; bar? [1=Yes]
;     ring:                               0L, $ ; ring? [1=Yes]
;     multiple:                           0L, $ ; multiple? [1=Yes]
;     sb25_b:                         -999.0, $ ; mean surface brightness within D25 (B-band)
;     ned_class:                       '...', $ ; AGN according to NED? [U/Y/N] (U=Unknown)
;     lit_class:                       '...', $ ; AGN according to the literature? [U/Y/N] (U=Unknown)
;     lit_class_ref:                   '...', $ ; LIT_CLASS reference
                              
      distance:                       -999.0, $ ; distance [Mpc]
      distance_err:                   -999.0, $ ; distance error [Mpc]
      distance_ref:                    '...', $ ; distance reference
      distance_texref:                 '...', $ ; distance reference (bibtex code)
      distance_method:                 '...', $ ; distance method
      distance_model:                 -999.0, $ ; model-based distance (Mould et al. 2000)
      ebv_mw:                         -999.0, $ ; Galactic E(B-V) (SFD)
                              
      d25_origin:                     '', $ ; origin
      d25_ref:                        '', $ ; reference
      d25_maj:                    -999.0, $ ; D25 major axis [arcmin]
      d25_min:                    -999.0, $ ; D25 minor axis [arcmin]
      optical_posangle_origin:        '', $ ; origin
      optical_posangle_ref:           '', $ ; reference
      optical_posangle:           -999.0, $ ; optical position angle (N-->E)
      optical_inclination:        -999.0, $ ; optical inclination angle
      twomass_maj:                -999.0, $ ; 2MASS total major axis diameter [arcmin]
      twomass_min:                -999.0, $ ; 2MASS total minor axis diameter [arcmin]
      twomass_k20_maj:            -999.0, $ ; 2MASS K20 major axis diameter [arcmin]
      twomass_k20_min:            -999.0, $ ; 2MASS K20 minor axis diameter [arcmin]
      twomass_diameter_ref:           '', $ ; reference
      twomass_posangle:           -999.0, $ ; infrared position angle (N-->E)
      twomass_inclination:        -999.0, $ ; infrared inclination angle
      posangle_origin:                '', $ ; origin
      posangle_ref:                   '', $ ; reference
      posangle:                   -999.0, $ ; "final" position angle (N-->E)
      inclination:                -999.0, $ ; "final" inclination angle

      rc3_ubv_ref:            replicate('...',3), $
      rc3_ubv:                   fltarr(3)-999.0, $
      rc3_ubv_err:               fltarr(3)-999.0, $
      rc3_ubv_absmag:            fltarr(3)-999.0, $
      rc3_ubv_absmag_err:        fltarr(3)-999.0, $

      twomass_jhk_ref:                     '...', $
      twomass_jhk:               fltarr(3)-999.0, $
      twomass_jhk_err:           fltarr(3)-999.0, $
      twomass_jhk_absmag:        fltarr(3)-999.0, $
      twomass_jhk_absmag_err:    fltarr(3)-999.0, $
      twomass_jhk20:             fltarr(3)-999.0, $
      twomass_jhk20_err:         fltarr(3)-999.0, $
      twomass_jhk20_absmag:      fltarr(3)-999.0, $
      twomass_jhk20_absmag_err:  fltarr(3)-999.0, $

      iras_ref:       replicate('...',4), $
      iras:              fltarr(4)-999.0, $
      iras_err:          fltarr(4)-999.0, $
      iras_absmag:       fltarr(4)-999.0, $
      iras_absmag_err:   fltarr(4)-999.0, $

      iras_12_predict:            -999.0, $
      iras_12_predict_err:        -999.0, $
      iras_25_predict:            -999.0, $
      iras_25_predict_err:        -999.0, $
      
      irac_ref:                    '...', $
      irac:              fltarr(4)-999.0, $
      irac_err:          fltarr(4)-999.0, $
      irac_absmag:       fltarr(4)-999.0, $
      irac_absmag_err:   fltarr(4)-999.0, $

      mips_ref:               '...', $
      mips:              fltarr(3)-999.0, $
      mips_err:          fltarr(3)-999.0, $
      mips_absmag:       fltarr(3)-999.0, $
      mips_absmag_err:   fltarr(3)-999.0, $

      galex_ref:              '...', $
      galex:             fltarr(2)-999.0, $
      galex_err:         fltarr(2)-999.0, $
      galex_absmag:      fltarr(2)-999.0, $
      galex_absmag_err:  fltarr(2)-999.0, $

      fir_flux_sm96:     fltarr(2)-999.0, $
      fir_lum_sm96:      fltarr(2)-999.0, $

      ir_flux:           fltarr(2)-999.0, $
      ir_lum:            fltarr(2)-999.0, $
      ir_dust_temp:               -999.0, $ ; dust temperature corresponding to the f(60)/f(100) ratio [K]
      ir_flux_sm96:      fltarr(2)-999.0, $
      ir_lum_sm96:       fltarr(2)-999.0, $
      ir_flux_dale01:    fltarr(2)-999.0, $
      ir_lum_dale01:     fltarr(2)-999.0, $
      ir_flux_dale02:    fltarr(2)-999.0, $
      ir_lum_dale02:     fltarr(2)-999.0}
      
;     U:                          -999.0, $ ; final U magnitude, corrected for emission lines
;     U_err:                      -999.0, $
;     U_lum:                      -999.0, $ ; U-band luminosity
;     U_lum_err:                  -999.0, $
;     M_U:                        -999.0, $ ; absolute U-band magnitude
;     M_U_err:                    -999.0, $
;     U_obs:                      -999.0, $ ; final U magnitude, un-corrected for emission lines
;     U_obs_err:                  -999.0, $
;     U_lum_obs:                  -999.0, $ ; U-band luminosity
;     U_lum_obs_err:              -999.0, $
;     M_U_obs:                    -999.0, $ ; absolute U-band magnitude
;     M_U_obs_err:                -999.0, $
;     B:                          -999.0, $ ; final B magnitude, corrected for emission lines
;     B_err:                      -999.0, $
;     B_lum:                      -999.0, $ ; B-band luminosity
;     B_lum_err:                  -999.0, $
;     M_B:                        -999.0, $ ; absolute B-band magnitude
;     M_B_err:                    -999.0, $
;     B_obs:                      -999.0, $ ; final B magnitude, un-corrected for emission lines
;     B_obs_err:                  -999.0, $
;     B_lum_obs:                  -999.0, $ ; B-band luminosity
;     B_lum_obs_err:              -999.0, $
;     M_B_obs:                    -999.0, $ ; absolute B-band magnitude
;     M_B_obs_err:                -999.0, $
;     V:                          -999.0, $ ; final V magnitude, corrected for emission lines
;     V_err:                      -999.0, $
;     V_lum:                      -999.0, $ ; V-band luminosity
;     V_lum_err:                  -999.0, $
;     M_V:                        -999.0, $ ; absolute V-band magnitude
;     M_V_err:                    -999.0, $
;     V_obs:                      -999.0, $ ; final V magnitude, un-corrected for emission lines
;     V_obs_err:                  -999.0, $
;     V_lum_obs:                  -999.0, $ ; V-band luminosity
;     V_lum_obs_err:              -999.0, $
;     M_V_obs:                    -999.0, $ ; absolute V-band magnitude
;     M_V_obs_err:                -999.0, $
;     R:                          -999.0, $ ; final R magnitude, corrected for emission lines
;     R_err:                      -999.0, $
;     R_lum:                      -999.0, $ ; R-band luminosity
;     R_lum_err:                  -999.0, $
;     M_R:                        -999.0, $ ; absolute R-band magnitude
;     M_R_err:                    -999.0, $
;     R_obs:                      -999.0, $ ; final R magnitude, un-corrected for emission lines
;     R_obs_err:                  -999.0, $
;     R_lum_obs:                  -999.0, $ ; R-band luminosity
;     R_lum_obs_err:              -999.0, $
;     M_R_obs:                    -999.0, $ ; absolute R-band magnitude
;     M_R_obs_err:                -999.0, $
;     I:                          -999.0, $ ; final I magnitude, corrected for emission lines
;     I_err:                      -999.0, $
;     I_lum:                      -999.0, $ ; I-band luminosity
;     I_lum_err:                  -999.0, $
;     M_I:                        -999.0, $ ; absolute I-band magnitude
;     M_I_err:                    -999.0, $
;     I_obs:                      -999.0, $ ; final I magnitude, un-corrected for emission lines
;     I_obs_err:                  -999.0, $
;     I_lum_obs:                  -999.0, $ ; I-band luminosity
;     I_lum_obs_err:              -999.0, $
;     M_I_obs:                    -999.0, $ ; absolute I-band magnitude
;     M_I_obs_err:                -999.0, $

;     rc3_u_origin:                   '', $ ; RC3 or LEDA?
;     rc3_U:                      -999.0, $ ; RC3/LEDA total apparent U magnitude and error
;     rc3_U_err:                  -999.0, $
;     rc3_U_lum:                  -999.0, $ ; U-band luminosity
;     rc3_U_lum_err:              -999.0, $
;     rc3_M_U:                    -999.0, $ ; absolute U-band magnitude
;     rc3_M_U_err:                -999.0, $
;     rc3_b_origin:                   '', $ ; RC3 or LEDA?
;     rc3_B:                      -999.0, $ ; RC3/LEDA total apparent B magnitude and error
;     rc3_B_err:                  -999.0, $
;     rc3_B_lum:                  -999.0, $ ; B-band luminosity
;     rc3_B_lum_err:              -999.0, $
;     rc3_M_B:                    -999.0, $ ; absolute B-band magnitude
;     rc3_M_B_err:                -999.0, $
;     rc3_v_origin:                   '', $ ; RC3 or LEDA?
;     rc3_V:                      -999.0, $ ; RC3/LEDA total apparent V magnitude and error
;     rc3_V_err:                  -999.0, $
;     rc3_V_lum:                  -999.0, $ ; V-band luminosity
;     rc3_V_lum_err:              -999.0, $
;     rc3_M_V:                    -999.0, $ ; absolute V-band magnitude
;     rc3_M_V_err:                -999.0, $

;     twomass_origin:                 '', $ ; 2MASS photometry origin [XSC or LGA]
;     twomass_J:                  -999.0, $ ; 2MASS J magnitude and error
;     twomass_J_err:              -999.0, $
;     twomass_J_lum:              -999.0, $ ; J-band luminosity
;     twomass_J_lum_err:          -999.0, $
;     twomass_M_J:                -999.0, $ ; absolute J-band magnitude
;     twomass_M_J_err:            -999.0, $
;     twomass_H:                  -999.0, $ ; 2MASS H magnitude and error
;     twomass_H_err:              -999.0, $
;     twomass_H_lum:              -999.0, $ ; H-band luminosity
;     twomass_H_lum_err:          -999.0, $
;     twomass_M_H:                -999.0, $ ; absolute H-band magnitude
;     twomass_M_H_err:            -999.0, $
;     twomass_Ks:                 -999.0, $ ; 2MASS Ks magnitude and error
;     twomass_Ks_err:             -999.0, $
;     twomass_Ks_lum:             -999.0, $ ; Ks-band luminosity
;     twomass_Ks_lum_err:         -999.0, $
;     twomass_M_Ks:               -999.0, $ ; absolute Ks-band magnitude
;     twomass_M_Ks_err:           -999.0, $
;     twomass_J20:                -999.0, $ ; 2MASS J20 magnitude and error
;     twomass_J20_err:            -999.0, $
;     twomass_J20_lum:            -999.0, $ ; J20-band luminosity
;     twomass_J20_lum_err:        -999.0, $
;     twomass_M_J20:              -999.0, $ ; absolute J20-band magnitude
;     twomass_M_J20_err:          -999.0, $
;     twomass_H20:                -999.0, $ ; 2MASS H20 magnitude and error
;     twomass_H20_err:            -999.0, $
;     twomass_H20_lum:            -999.0, $ ; H20-band luminosity
;     twomass_H20_lum_err:        -999.0, $
;     twomass_M_H20:              -999.0, $ ; absolute H20-band magnitude
;     twomass_M_H20_err:          -999.0, $
;     twomass_Ks20:               -999.0, $ ; 2MASS Ks20 magnitude and error
;     twomass_Ks20_err:           -999.0, $
;     twomass_Ks20_lum:           -999.0, $ ; Ks20-band luminosity
;     twomass_Ks20_lum_err:       -999.0, $
;     twomass_M_Ks20:             -999.0, $ ; absolute Ks20-band magnitude
;     twomass_M_Ks20_err:         -999.0, $

;     fuv_FLUX:                   -999.0, $ ; FUV flux from Bell (2003) [erg/s/cm2]
;     fuv_FLUX_ERR:               -999.0, $
;     fuv_LUM:                    -999.0, $ ; FUV luminosity [erg/s]
;     fuv_LUM_ERR:                -999.0, $
;     fuv_REF:                       ' ', $
;     radio_1_4GHZ_FLUX:          -999.0, $ ; 1.4 GHz flux from Bell (2003) [erg/s/cm2]
;     radio_1_4GHZ_FLUX_ERR:      -999.0, $
;     radio_1_4GHZ_LUM:           -999.0, $ ; 1.4 GHz luminosity [erg/s]
;     radio_1_4GHZ_LUM_ERR:       -999.0, $
;     W20:                        -999.0, $ ; 21 cm line width at 20% of peak (km/s)
;     W20_err:                    -999.0, $
;     mass_HI_leda:               -999.0, $ ; LEDA HI mass corrected for self-absorption (log M_sun)
;     mass_HI_leda_err:           -999.0, $
;     mass_HI_lit:                -999.0, $ ; Literature HI mass (log M_sun)
;     mass_HI_lit_err:            -999.0, $
;     mass_HI_lit_ref:               ' ', $
;     mass_CO_lit:                -999.0, $ ; Literature CO mass (log M_sun)
;     mass_CO_lit_err:            -999.0, $
;     mass_CO_lit_ref:               ' ', $
;     eta_gas:                    -999.0, $ ; Equation (4) in McGaugh & de Block (1997)
;     mass_gas:                   -999.0, $ ; total gas mass
;     mass_gas_err:               -999.0, $

;     iras_12:                    -999.0, $ ; IRAS 12 micron flux and error [Jy] (upper limits are negative)
;     iras_12_err:                -999.0, $
;     iras_12_ref:                   ' ', $
;     iras_12_predict:            -999.0, $ ; predicted IRAS 12 micron flux and error [Jy]
;     iras_12_predict_err:        -999.0, $
;     iras_25:                    -999.0, $ ; IRAS 25 micron flux and error [Jy]
;     iras_25_err:                -999.0, $
;     iras_25_ref:                   ' ', $
;     iras_25_predict:            -999.0, $ ; predicted IRAS 25 micron flux and error [Jy]
;     iras_25_predict_err:        -999.0, $
;     iras_60:                    -999.0, $ ; IRAS 60 micron flux and error [Jy]
;     iras_60_err:                -999.0, $
;     iras_60_ref:                   ' ', $
;     iras_100:                   -999.0, $ ; IRAS 100 micron flux and error [Jy]
;     iras_100_err:               -999.0, $
;     iras_100_ref:                  ' ', $

;     fir_flux:                   -999.0, $ ; FIR flux (Helou et al 1988) [erg/s/cm2]
;     fir_flux_err:               -999.0, $
;     fir_lum:                    -999.0, $ ; FIR luminosity [erg/s]
;     fir_lum_err:                -999.0, $ 
 ;    ir_flux:                    -999.0, $ ; total IR flux (direct integration) [erg/s/cm2]
;     ir_flux_err:                -999.0, $
;     ir_dust_temp:               -999.0, $ ; dust temperature corresponding to the f(60)/f(100) ratio [K]
;     ir_flux_dale01:             -999.0, $ ; total IR flux (Dale et al. 2001) [erg/s/cm2]
;     ir_flux_dale01_err:         -999.0, $
;     ir_flux_dale02:             -999.0, $ ; total IR flux (Dale & Helou 2002) [erg/s/cm2]
;     ir_flux_dale02_err:         -999.0, $
;     ir_flux_sm96:               -999.0, $ ; total IR flux (Sanders & Mirabel 1996) [erg/s/cm2]
;     ir_flux_sm96_err:           -999.0, $
;     ir_lum:                     -999.0, $ ; total IR luminosity [erg/s]
;     ir_lum_err:                 -999.0, $
;     ir_lum_dale02:              -999.0, $ ; total IR luminosity [erg/s]
;     ir_lum_dale02_err:          -999.0, $
;     ir_lum_sm96:                -999.0, $ ; total IR luminosity [erg/s]
;     ir_lum_sm96_err:            -999.0, $
;     ir_lum_dale01:              -999.0, $ ; total IR luminosity [erg/s]
;     ir_lum_dale01_err:          -999.0}
;     sfr_ir:                     -999.0, $ ; IR star-formation rate (Kennicutt 1998)
;     sfr_ir_err:                 -999.0, $
;     sfr_ir_dale01:              -999.0, $ ; IR star-formation rate (Kennicutt 1998)
;     sfr_ir_dale01_err:          -999.0, $
;     L_FIR_L_B:                  -999.0, $ ; FIR-to-optical luminosity ratio and error
;     L_FIR_L_B_err:              -999.0, $
;     L_IR_L_B:                   -999.0, $ ; IR-to-optical luminosity ratio and error
;     L_IR_L_B_err:               -999.0, $

;     A_U_incl:                   -999.0, $ ; extinction correction to face-on inclination
;     A_B_incl:                   -999.0, $
;     A_V_incl:                   -999.0}
;     lit_log12oh:                -999.0, $ ; literature 12+log(O/H) (electron temperature) and error
;     lit_log12oh_err:            -999.0, $
;     lit_log12oh_ref:               ' '}
    table = replicate(table,ngalaxy)

return, table
end
