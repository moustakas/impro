pro parse_storchi95, out
; jm06dec13nyu 

    outpath = hiiregions_path()

    data = rsex(outpath+'95storchi/raw_storchi95.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_one_hiiregion_structure(nobject)

    out.hii_galaxy = strcompress(data.galaxy,/remove)
    out.hii_region = 'H1'

; correct the Balmer line for stellar absorption
    niceprint, data.galaxy, data.hb_abs*data.hb_cont/data.h_beta
    data.h_beta = data.h_beta + data.hb_abs*data.hb_cont
    
; reddening-correct the data
    line = {$
      linename:         ['OII_3727','H_BETA','OIII_5007','H_ALPHA','NII_6584','SII_6716','SII_6731'], $
      oii_3727_wave:    3727.42, $
      h_beta_wave:      4861.33, $
      oiii_5007_wave:   5006.84, $
      h_alpha_wave:     6562.80, $
      nii_6584_wave:    6583.46, $
      sii_6716_wave:    6716.14, $
      sii_6731_wave:    6730.81, $
      oii_3727:      [0.0,-2.0], $
      h_beta:        [0.0,-2.0], $
      oiii_5007:     [0.0,-2.0], $
      h_alpha:       [0.0,-2.0], $
      nii_6584:      [0.0,-2.0], $
      sii_6716:      [0.0,-2.0], $
      sii_6731:      [0.0,-2.0]}
    line = replicate(line,nobject)

    k_oii_3727  = k_lambda(line[0].oii_3727_wave,/odonnell)
    k_h_beta    = k_lambda(line[0].h_beta_wave,/odonnell)
    k_oiii_5007 = k_lambda(line[0].oiii_5007_wave,/odonnell)
    k_h_alpha   = k_lambda(line[0].h_alpha_wave,/odonnell)
    k_nii_6584  = k_lambda(line[0].nii_6584_wave,/odonnell)
    k_sii_6716  = k_lambda(line[0].sii_6716_wave,/odonnell)
    k_sii_6731  = k_lambda(line[0].sii_6731_wave,/odonnell)
    
    line.oii_3727[0]  = data.oii_3727*1D-15
    line.h_beta[0]    = data.h_beta*1D-15
    line.oiii_5007[0] = data.oiii_5007*1D-15
    line.h_alpha[0]   = data.h_alpha*1D-15
    line.nii_6584[0]  = data.nii_6584*1D-15
    line.sii_6716[0]  = data.sii_6716*1D-15
    line.sii_6731[0]  = data.sii_6731*1D-15

    line.oii_3727[1]  = data.oii_3727_err*1D-15
    line.h_beta[1]    = data.h_beta_err*1D-15
    line.oiii_5007[1] = data.oiii_5007_err*1D-15
    line.h_alpha[1]   = data.h_alpha_err*1D-15
    line.nii_6584[1]  = data.nii_6584_err*1D-15
    line.sii_6716[1]  = data.sii_6716_err*1D-15
    line.sii_6731[1]  = data.sii_6731_err*1D-15

; correct for reddening    
    
    line_nodust = iunred_linedust(line,/silent,/nopropagate)
    niceprint, data.galaxy, line_nodust.hahb, line_nodust.hahb_err, $
      line_nodust.ebv_hahb, line_nodust.ebv_hahb_err
    

; construct the final line ratios, relative to H-beta
    
    out.oii_3727 = line_nodust.oii_3727[0]/line_nodust.h_beta[0]
    out.oii_3727_err = im_compute_error(line_nodust.oii_3727[0],line_nodust.oii_3727[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

    out.oiii_5007 = line_nodust.oiii_5007[0]/line_nodust.h_beta[0]
    out.oiii_5007_err = im_compute_error(line_nodust.oiii_5007[0],line_nodust.oiii_5007[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

    out.ha = line_nodust.h_alpha[0]/line_nodust.h_beta[0]
    out.ha_err = im_compute_error(line_nodust.h_alpha[0],line_nodust.h_alpha[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

    out.nii_6584 = line_nodust.nii_6584[0]/line_nodust.h_beta[0]
    out.nii_6584_err = im_compute_error(line_nodust.nii_6584[0],line_nodust.nii_6584[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

    out.sii_6716 = line_nodust.sii_6716[0]/line_nodust.h_beta[0]
    out.sii_6716_err = im_compute_error(line_nodust.sii_6716[0],line_nodust.sii_6716[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

    out.sii_6731 = line_nodust.sii_6731[0]/line_nodust.h_beta[0]
    out.sii_6731_err = im_compute_error(line_nodust.sii_6731[0],line_nodust.sii_6731[1],$
      line_nodust.h_beta[0],line_nodust.h_beta[1],/quotient)

; write out

    filename = '1995_storchi.sex'
    reference = 'Storchi-Bergmann et al. 1995, ApJS, 98, 103'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_STORCHI95 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    printf, lun, '## Only a subset of the available data have been tabulated.'
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
