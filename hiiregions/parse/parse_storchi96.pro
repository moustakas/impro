pro parse_storchi96, out
; jm06jan18uofa - parse Table 3 from Hunter & Hoffman 1999

    outpath = hiiregions_path()

    data = rsex(outpath+'96storchi/raw_storchi96.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.hii_galaxy,/remove)
    out.hii_region = strcompress(data.hii_region,/remove)

    nw = where(data.radius gt 0.0,comp=se)
    out[nw].raoffset = -data[nw].radius*cos((data[nw].pa-90.0)*!dtor)
    out[nw].deoffset = +data[nw].radius*sin((data[nw].pa-90.0)*!dtor)
    out[se].raoffset = -data[se].radius*cos((data[se].pa-90.0)*!dtor) > 0.0
    out[se].deoffset = +data[se].radius*sin((data[se].pa-90.0)*!dtor)

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

    k_h_beta    = k_lambda(line[0].h_beta_wave,/odonnell)
    k_h_alpha   = k_lambda(line[0].h_alpha_wave,/odonnell)
    k_oiii_5007 = k_lambda(line[0].oiii_5007_wave,/odonnell)
    k_nii_6584  = k_lambda(line[0].nii_6584_wave,/odonnell) 
    k_sii_6716  = k_lambda(line[0].sii_6716_wave,/odonnell) 
    k_sii_6731  = k_lambda(line[0].sii_6731_wave,/odonnell) 
    
    line.oii_3727[0] = data.oii_3727*1D-15
    lo = where(line.oii_3727[0] lt 1D-15,comp=hi)
    line[lo].oii_3727[1] = line[lo].oii_3727[0]*0.2
    line[hi].oii_3727[1] = line[hi].oii_3727[0]*0.1
    
    line.h_alpha[0] = data.h_alpha*1D-15
    lo = where(line.h_alpha[0] lt 1D-15,comp=hi)
    line[hi].h_alpha[1] = line[hi].h_alpha[0]*0.1

; predict the observed, absorption-corrected H-beta flux based on the
; reddening

    line.h_beta[0] = line.h_alpha[0]*10^(-0.4*(k_h_beta-k_h_alpha)*data.ebv)/2.86
    line.h_beta[1] = line.h_alpha[1]*10^(-0.4*(k_h_beta-k_h_alpha)*data.ebv)/2.86
    
;   line.h_beta[0] = data.h_beta*1D-15
;   lo = where(line.h_beta[0] lt 1D-15,comp=hi)
;   line[lo].h_beta[1] = line[lo].h_beta[0]*0.2
;   line[hi].h_beta[1] = line[hi].h_beta[0]*0.1

    line.oiii_5007[0] = data.oiii_5007*1D-15
    lo = where(line.oiii_5007[0] lt 1D-15,comp=hi)
    line[lo].oiii_5007[1] = line[lo].oiii_5007[0]*0.2
    line[hi].oiii_5007[1] = line[hi].oiii_5007[0]*0.1

    line.nii_6584[0] = data.nii_6584*1D-15
    lo = where(line.nii_6584[0] lt 1D-15,comp=hi)
    line[hi].nii_6584[1] = line[hi].nii_6584[0]*0.1

    line.sii_6716[0] = data.sii_6716*1D-15
    lo = where(line.sii_6716[0] lt 1D-15,comp=hi)
    line[lo].sii_6716[1] = line[lo].sii_6716[0]*0.2
    line[hi].sii_6716[1] = line[hi].sii_6716[0]*0.1

    line.sii_6731[0] = data.sii_6731*1D-15
    lo = where(line.sii_6731[0] lt 1D-15,comp=hi)
    line[lo].sii_6731[1] = line[lo].sii_6731[0]*0.2
    line[hi].sii_6731[1] = line[hi].sii_6731[0]*0.1

    line_nodust = iunred_linedust(line,/silent,/nopropagate)
;   niceprint, line_nodust.ebv_hahb, data.ebv
    
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

    filename = '1996_storchi.sex'
    reference = 'Storchi-Bergmann et al. 1996, ApJ, 460, 252'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_STORCHI96 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    printf, lun, '## Data on NGC1672 and NGC5248 not tabulated.'
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
