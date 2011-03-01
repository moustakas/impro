pro parse_kong99, out
; jm06dec13nyu 

    outpath = hiiregions_path()

    data = rsex(outpath+'99kong/raw_kong99.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.galaxy,/remove)
    out.hii_region = 'H1'

; reddening-correct the data

    k_oii_3727  = k_lambda(3727.42,/odonnell)
    k_h_beta    = k_lambda(4861.33,/odonnell)
    k_oiii_5007 = k_lambda(5006.84,/odonnell)
    k_h_alpha   = k_lambda(6562.80,/odonnell)
    k_nii_6584  = k_lambda(6583.46,/odonnell)
    
    ebv = get_ebv(1.0/data.h_beta)

    h_beta = data.h_beta*10.0^(0.4*ebv*(k_h_beta-k_h_alpha))
    
    out.ha = 1.0/h_beta
    out.ha_err = 0.05*out.ha

    out.oii_3727     = data.oii_3727*10.0^(0.4*ebv*(k_oii_3727-k_h_alpha))/h_beta ; [OII]/Hb
    out.oii_3727_err = 0.05*out.oii_3727

    out.oiii_5007     = data.oiii_5007*10.0^(0.4*ebv*(k_oiii_5007-k_h_alpha))/h_beta ; [OIII]/Hb
    out.oiii_5007_err = 0.05*out.oiii_5007
    
    out.nii_6584     = data.nii_6584*10.0^(0.4*ebv*(k_nii_6584-k_h_alpha))/h_beta ; [NII]/Hb
    out.nii_6584_err = 0.05*out.nii_6584

; write out

    filename = '1999_kong.sex'
    reference = 'Kong & Cheng 1999, A&A, 351, 477'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_KONG99 '+im_today()+'.'
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
