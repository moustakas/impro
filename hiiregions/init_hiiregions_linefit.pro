function init_hiiregions_linefit, data
; generate a line flux structure equivalent to the output from
; ISPECLINEFIT() that can be used in IM_ABUNDANCE() and
; WRITE_HII_REGIONS(); all the line fluxes are assumed to be reddening
; corrected and normalized to H-beta ( = 1.0)

; jm06sep28nyu - excised from WRITE_HII_REGIONS    
    
    nobject = n_elements(data)
    
    linefit = {$
      hii_region:           '', $
      texref:               '', $
      linename:     strarr(15), $ ; this structure tricks IM_ABUNDANCE()
      ebv_hahb:     0.0,        $
      ebv_hahb_err: 0.1,        $ 
      oii_3727:     [0.0,-2.0], $
      oiii_4363:    [0.0,-2.0], $
      h_beta:       [0.0,-2.0], $
      oiii_4959:    [0.0,-2.0], $
      oiii_5007:    [0.0,-2.0], $
      nii_5755:     [0.0,-2.0], $
      siii_6312:    [0.0,-2.0], $
      nii_6548:     [0.0,-2.0], $
      h_alpha:      [0.0,-2.0], $
      nii_6584:     [0.0,-2.0], $
      sii_6716:     [0.0,-2.0], $
      sii_6731:     [0.0,-2.0], $
      oii_7325:     [0.0,-2.0], $
      siii_9069:    [0.0,-2.0], $
      siii_9532:    [0.0,-2.0]}
    linefit = replicate(linefit,nobject)
    linefit.linename = ['OII_3727','OIII_4363','H_BETA','OIII_4959',$
      'OIII_5007','NII_5755','SIII_6312','NII_6548','H_ALPHA','NII_6584',$
      'SII_6716','SII_6731','OII_7325','SIII_9069','SIII_9532']

    linefit.hii_region = data.hii_region
    linefit.texref = data.texref
    
    linefit.h_beta = [1.0,0.001]
;   linefit.h_alpha = [return_tbalmer(),0.001] ; T = 10000 K

    g = where(data.h_alpha_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].h_alpha = transpose([ [data[g].h_alpha_h_beta], [data[g].h_alpha_h_beta_err] ])

    g = where(data.oii_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].oii_3727 = transpose([ [data[g].oii_h_beta], [data[g].oii_h_beta_err] ])

    g = where(data.oiii_4363_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].oiii_4363 = transpose([ [data[g].oiii_4363_h_beta], [data[g].oiii_4363_h_beta_err] ])

    g = where(data.oiii_4959_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].oiii_4959 = transpose([ [data[g].oiii_4959_h_beta], [data[g].oiii_4959_h_beta_err] ])

    g = where(data.oiii_5007_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].oiii_5007 = transpose([ [data[g].oiii_5007_h_beta], [data[g].oiii_5007_h_beta_err] ])

    g = where(data.nii_5755_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].nii_5755 = transpose([ [data[g].nii_5755_h_beta], [data[g].nii_5755_h_beta_err] ])

    g = where(data.siii_6312_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].siii_6312 = transpose([ [data[g].siii_6312_h_beta], [data[g].siii_6312_h_beta_err] ])

    g = where(data.nii_6548_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].nii_6548 = transpose([ [data[g].nii_6548_h_beta], [data[g].nii_6548_h_beta_err] ])

    g = where(data.nii_6584_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].nii_6584 = transpose([ [data[g].nii_6584_h_beta], [data[g].nii_6584_h_beta_err] ])

    g = where(data.sii_6716_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].sii_6716 = transpose([ [data[g].sii_6716_h_beta], [data[g].sii_6716_h_beta_err] ])

    g = where(data.sii_6731_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].sii_6731 = transpose([ [data[g].sii_6731_h_beta], [data[g].sii_6731_h_beta_err] ])

    g = where(data.oii_7325_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].oii_7325 = transpose([ [data[g].oii_7325_h_beta], [data[g].oii_7325_h_beta_err] ])

    g = where(data.siii_9069_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].siii_9069 = transpose([ [data[g].siii_9069_h_beta], [data[g].siii_9069_h_beta_err] ])

    g = where(data.siii_9532_h_beta gt -900.0,ng)
    if (ng ne 0L) then linefit[g].siii_9532 = transpose([ [data[g].siii_9532_h_beta], [data[g].siii_9532_h_beta_err] ])

return, linefit
end

