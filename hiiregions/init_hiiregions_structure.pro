function init_hiiregions_structure, nhii
; initialize the HII region data structure
; jm07nov26nyu - excised from WRITE_HII_REGIONS
    
    hii = {$
      hii_galaxy:                                ' ', $
      ned_galaxy:                                ' ', $
      hii_region:                                ' ', $
                                             
      galaxy_ra:                                 ' ', $
      galaxy_dec:                                ' ', $
      galaxy_rc3_r25:                         -999.0, $ ; RC3 major axis radius [arcmin]
      galaxy_rc3_pa:                          -999.0, $ ; RC3 galaxy position angle
      galaxy_rc3_incl:                        -999.0, $ ; RC3 galaxy inclination
      galaxy_twomass_k20:                     -999.0, $ ; 2MASS major axis radius [arcmin]
      galaxy_twomass_pa:                      -999.0, $ ; 2MASS galaxy position angle
      galaxy_twomass_incl:                    -999.0, $ ; 2MASS galaxy inclination

      hii_raoffset:                           -999.0, $
      hii_deoffset:                           -999.0, $
      hii_rc3_radius:                         -999.0, $ ; de-projected galactocentric radius [arcmin]
      hii_rc3_rr25:                           -999.0, $ ; de-projected galactocentric R/R25
      hii_rc3_phi:                            -999.0, $ ; HII-region position angle
      hii_twomass_radius:                     -999.0, $ ; de-projected galactocentric radius [arcmin]
      hii_twomass_rr25:                       -999.0, $ ; de-projected galactocentric R/R25
      hii_twomass_phi:                        -999.0, $ ; HII-region position angle
      hii_lit_radius:                         -999.0, $ ; literature de-projected galactocentric radius [arcmin]
      hii_lit_rr25:                           -999.0, $ ; literature de-projected galactocentric R/R25

      texref:                                     '', $
      reference:                                  '', $

      oii_h_alpha:                           -999.0, $ ; 3727
      oii_h_alpha_err:                       -999.0, $
      oii_h_beta:                            -999.0, $
      oii_h_beta_err:                        -999.0, $
                                             
      oii_7325_h_alpha:                      -999.0, $ ; 7325
      oii_7325_h_alpha_err:                  -999.0, $
      oii_7325_h_beta:                       -999.0, $
      oii_7325_h_beta_err:                   -999.0, $
                                             
      oiii_4363_h_alpha:                     -999.0, $ ; 4363
      oiii_4363_h_alpha_err:                 -999.0, $
      oiii_4363_h_beta:                      -999.0, $
      oiii_4363_h_beta_err:                  -999.0, $
                                             
      oiii_4959_h_alpha:                     -999.0, $ ; 4959
      oiii_4959_h_alpha_err:                 -999.0, $
      oiii_4959_h_beta:                      -999.0, $
      oiii_4959_h_beta_err:                  -999.0, $
                                             
      oiii_5007_h_alpha:                     -999.0, $ ; 5007
      oiii_5007_h_alpha_err:                 -999.0, $
      oiii_5007_h_beta:                      -999.0, $
      oiii_5007_h_beta_err:                  -999.0, $
                                             
      oiii_h_alpha:                          -999.0, $ ; 4959+5007
      oiii_h_alpha_err:                      -999.0, $
      oiii_h_beta:                           -999.0, $
      oiii_h_beta_err:                       -999.0, $
                                             
      nii_5755_h_alpha:                      -999.0, $ ; 5755
      nii_5755_h_alpha_err:                  -999.0, $
      nii_5755_h_beta:                       -999.0, $
      nii_5755_h_beta_err:                   -999.0, $

      nii_6548_h_alpha:                      -999.0, $ ; 6548
      nii_6548_h_alpha_err:                  -999.0, $
      nii_6548_h_beta:                       -999.0, $
      nii_6548_h_beta_err:                   -999.0, $

      h_alpha_h_beta:                        -999.0, $
      h_alpha_h_beta_err:                    -999.0, $
      
      nii_6584_h_alpha:                      -999.0, $ ; 6584
      nii_6584_h_alpha_err:                  -999.0, $
      nii_6584_h_beta:                       -999.0, $
      nii_6584_h_beta_err:                   -999.0, $
                                             
      nii_h_alpha:                           -999.0, $ ; 6548+6584
      nii_h_alpha_err:                       -999.0, $
      nii_h_beta:                            -999.0, $
      nii_h_beta_err:                        -999.0, $
                                             
      sii_6716_h_alpha:                      -999.0, $ ; 6716
      sii_6716_h_alpha_err:                  -999.0, $
      sii_6716_h_beta:                       -999.0, $
      sii_6716_h_beta_err:                   -999.0, $
                                             
      sii_6731_h_alpha:                      -999.0, $ ; 6731
      sii_6731_h_alpha_err:                  -999.0, $
      sii_6731_h_beta:                       -999.0, $
      sii_6731_h_beta_err:                   -999.0, $
                                             
      siii_6312_h_alpha:                     -999.0, $ ; 6312
      siii_6312_h_alpha_err:                 -999.0, $
      siii_6312_h_beta:                      -999.0, $
      siii_6312_h_beta_err:                  -999.0, $
                                             
      siii_9069_h_alpha:                     -999.0, $ ; 9069
      siii_9069_h_alpha_err:                 -999.0, $
      siii_9069_h_beta:                      -999.0, $
      siii_9069_h_beta_err:                  -999.0, $
                                             
      siii_9532_h_alpha:                     -999.0, $ ; 9532
      siii_9532_h_alpha_err:                 -999.0, $
      siii_9532_h_beta:                      -999.0, $
      siii_9532_h_beta_err:                  -999.0, $
                                             
      siii_h_alpha:                          -999.0, $ ; 9069+9532
      siii_h_alpha_err:                      -999.0, $
      siii_h_beta:                           -999.0, $
      siii_h_beta_err:                       -999.0, $
                                             
      sii_h_alpha:                           -999.0, $ ; 6716+6731
      sii_h_alpha_err:                       -999.0, $
      sii_h_beta:                            -999.0, $
      sii_h_beta_err:                        -999.0, $
                                             
      nii_6584_oii:                          -999.0, $ ; 6584/3727
      nii_6584_oii_err:                      -999.0, $
      nii_6584_oiii_5007:                    -999.0, $ ; 6584/5007
      nii_6584_oiii_5007_err:                -999.0, $
      nii_6584_sii:                          -999.0, $ ; 6584/(6716+6731)
      nii_6584_sii_err:                      -999.0, $
                                             
      oii_nii_6584:                          -999.0, $ ; 3727/6584
      oii_nii_6584_err:                      -999.0, $
      oii_oiii_5007:                         -999.0, $ ; 3727/5007
      oii_oiii_5007_err:                     -999.0, $
      oii_sii:                               -999.0, $ ; 3727/(6716+6731)
      oii_sii_err:                           -999.0, $
      oii_h_beta_oii_sii:                    -999.0, $ ; [3727/4861]*[3727/(6716+6731)]
      oii_h_beta_oii_sii_err:                -999.0, $ 
                                             
      oiii_5007_nii_6584:                    -999.0, $ ; 5007/6584
      oiii_5007_nii_6584_err:                -999.0, $ 
      oiii_5007_oii:                         -999.0, $ ; 5007/3727
      oiii_5007_oii_err:                     -999.0, $ 
      oiii_5007_sii:                         -999.0, $ ; 5007/(6716+6731)
      oiii_5007_sii_err:                     -999.0, $ 
      oiii_5007_h_beta_nii_6584_h_alpha:     -999.0, $ ; Pettini & Pagel 2004
      oiii_5007_h_beta_nii_6584_h_alpha_err: -999.0, $
                             
;     o32:                                   -999.0, $ ; (4959+5007)/3727
;     o32_err:                               -999.0, $
;     r23:                                   -999.0, $ ; (3727+4959+5007)/H-beta
;     r23_err:                               -999.0, $
;     s23:                                   -999.0, $ ; (6716+6731+9069+9532)/H-beta
;     s23_err:                               -999.0, $
                                             
      lit_t4363:                             -999.0, $ ; literature T(4363)
      lit_t4363_err:                         -999.0, $
      lit_t5755:                             -999.0, $ ; literature T(5755)
      lit_t5755_err:                         -999.0, $
      lit_t6312:                             -999.0, $ ; literature T(6312)
      lit_t6312_err:                         -999.0, $
      lit_t7325:                             -999.0, $ ; literature T(7325)
      lit_t7325_err:                         -999.0, $

;     lit_toiii:                             -999.0, $ ; literature T(O++)
;     lit_toiii_err:                         -999.0, $
;     lit_toii:                              -999.0, $ ; literature T(O+)
;     lit_toii_err:                          -999.0, $
;     lit_tnii:                              -999.0, $ ; literature T(N+)
;     lit_tnii_err:                          -999.0, $
;     lit_tsiii:                             -999.0, $ ; literature T(S++)
;     lit_tsiii_err:                         -999.0, $

      lit_log12oh_te:                        -999.0, $ ; literature electron-temperature abundance
      lit_log12oh_te_err:                    -999.0, $
                                             
      ewhb:                                  -999.0, $ ; EW(Hb) [Angstrom]
      fha:                                   -999.0}   ; f(Ha) [erg/s/cm2]
    hii = replicate(hii,nhii)

return, hii
end
