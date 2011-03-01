function init_one_hiiregion_structure, nobject
; jm05may12uofa
    
    if (n_elements(nobject) eq 0L) then nobject = 1L
    
    out = {$
      hii_galaxy:             '', $
      hii_region:             '', $
      raoffset:           -999.0, $
      deoffset:           -999.0, $
      radius:             -999.0, $
      rr25:               -999.0, $
      oii_3727:           -999.0, $
      oii_3727_err:       -999.0, $
      oiii_4363:          -999.0, $
      oiii_4363_err:      -999.0, $
      oiii_5007:          -999.0, $
      oiii_5007_err:      -999.0, $
      nii_5755:           -999.0, $
      nii_5755_err:       -999.0, $
      siii_6312:          -999.0, $
      siii_6312_err:      -999.0, $
      ha:                 -999.0, $
      ha_err:             -999.0, $
      nii_6584:           -999.0, $
      nii_6584_err:       -999.0, $
      sii_6716:           -999.0, $
      sii_6716_err:       -999.0, $
      sii_6731:           -999.0, $
      sii_6731_err:       -999.0, $
      oii_7325:           -999.0, $ ; = 7320,7330
      oii_7325_err:       -999.0, $
      siii_9069:          -999.0, $
      siii_9069_err:      -999.0, $
      siii_9532:          -999.0, $
      siii_9532_err:      -999.0, $
      siii_9069_9532:     -999.0, $
      siii_9069_9532_err: -999.0, $
      ewhb:               -999.0, $
      fha:                -999.0, $
      t4363:              -999.0, $
      t4363_err:          -999.0, $
      t5755:              -999.0, $
      t5755_err:          -999.0, $
      t6312:              -999.0, $
      t6312_err:          -999.0, $
      t7325:              -999.0, $
      t7325_err:          -999.0, $
;     toiii:              -999.0, $
;     toiii_err:          -999.0, $
;     toii:               -999.0, $
;     toii_err:           -999.0, $
;;    tnii:               -999.0, $ ; don't use this; always assume T(O+)=T(N+)
;;    tnii_err:           -999.0, $
;     tsiii:              -999.0, $
;     tsiii_err:          -999.0, $
      log12oh:            -999.0, $
      log12oh_err:        -999.0}
    out = replicate(out,nobject)

return, out
end
