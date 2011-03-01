pro parse_skillman94
; jm05may12uofa - parse Tables 2, 3, and 4 from Skillman et al. 1994
; jm07nov05nyu - separate out the three instruments into three
;                independent measurements

    outpath = hiiregions_path()

    nobject = 3L
    out = init_hii_region_structure(nobject)

    out.hii_galaxy = 'UGC04483'
    out.hii_region = ['INT','McD','WHT']

; tabulate the fluxes and errors    
    
    out.oii_3727 = [(0.337+0.498),0.944,0.837]
    out.oii_3727_err = [sqrt(0.023^2+0.034^2),0.120,0.062]

    out.oiii_4363 = [0.066,0.063,0.071]
    out.oiii_4363_err = [0.004,0.009,0.007]

    out.oiii_5007 = [2.690,2.530,3.070]
    out.oiii_5007_err = [0.07,0.08,0.12]

    out.ha = [2.780,2.780,2.780]
    out.ha_err = [0.11,0.06,0.075]

    out.nii_6584 = [-999.0,0.0317,0.029]
    out.nii_6584_err = [-999.0,0.0008,0.001]

    out.sii_6716 = [0.074,0.0726,0.069]
    out.sii_6716_err = [0.012,0.0017,0.002]

    out.sii_6731 = [0.051,0.0517,0.049]
    out.sii_6731_err = [0.008,0.0012,0.002]

    out.t4363 = [1.67,1.69,1.63]
    out.t4363_err = [0.06,0.13,0.08]

    out.log12oh = [7.50,7.48,7.57]
    out.log12oh_err = [0.04,0.06,0.05]

    out.ewhb = [170.0,165.0,144.0]
    
; write out

    filename = '1994_skillman.sex'
    reference = 'Skillman et al. 1994, ApJ, 431, 172'
    comments = 'The observations from each instrument were taken to be independent.'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_SKILLMAN94 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '+comments
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments

return
end
    

;;    nobject = 1L
;;    out = init_hii_region_structure(nobject)
;;
;;    out.hii_galaxy = 'UGC04483'
;;    out.hii_region = 'H1'
;;
;;; tabulate the fluxes and errors    
;;    
;;    oii = [(0.337+0.498),0.944,0.837]
;;    oii_err = [sqrt(0.023^2+0.034^2),0.120,0.062]
;;
;;    oiii_4363 = [0.066,0.063,0.071]
;;    oiii_4363_err = [0.004,0.009,0.007]
;;
;;    oiii = [2.690,2.530,3.070]
;;    oiii_err = [0.07,0.08,0.12]
;;
;;    ha = [2.780,2.780,2.780]
;;    ha_err = [0.11,0.06,0.075]
;;
;;    nii = [0.0317,0.029]
;;    nii_err = [0.0008,0.001]
;;
;;    sii_6716 = [0.074,0.0726,0.069]
;;    sii_6716_err = [0.012,0.0017,0.002]
;;
;;    sii_6731 = [0.051,0.0517,0.049]
;;    sii_6731_err = [0.008,0.0012,0.002]
;;
;;    toiii = [1.67,1.69,1.63]
;;    toiii_err = [0.06,0.13,0.08]
;;
;;    oh = [7.50,7.48,7.57]
;;    oh_err = [0.04,0.06,0.05]
;;
;;; compute the weighted mean fluxes
;;    
;;    out.oii_3727 = total(oii/oii_err^2)/total(1.0/oii_err^2)
;;    out.oii_3727_err = 1.0/sqrt(total(1.0/oii_err^2))
;;
;;    out.oiii_4363 = total(oiii_4363/oiii_4363_err^2)/total(1.0/oiii_4363_err^2)
;;    out.oiii_4363_err = 1.0/sqrt(total(1.0/oiii_4363_err^2))
;;
;;    out.oiii_5007 = total(oiii/oiii_err^2)/total(1.0/oiii_err^2)
;;    out.oiii_5007_err = 1.0/sqrt(total(1.0/oiii_err^2))
;;
;;    out.ha = total(ha/ha_err^2)/total(1.0/ha_err^2)
;;    out.ha_err = 1.0/sqrt(total(1.0/ha_err^2))
;;
;;    out.nii_6584 = total(nii/nii_err^2)/total(1.0/nii_err^2)
;;    out.nii_6584_err = 1.0/sqrt(total(1.0/nii_err^2))
;;
;;    out.sii_6716 = total(sii_6716/sii_6716_err^2)/total(1.0/sii_6716_err^2)
;;    out.sii_6716_err = 1.0/sqrt(total(1.0/sii_6716_err^2))
;;
;;    out.sii_6731 = total(sii_6731/sii_6731_err^2)/total(1.0/sii_6731_err^2)
;;    out.sii_6731_err = 1.0/sqrt(total(1.0/sii_6731_err^2))
;;
;;    out.toiii = total(toiii/toiii_err^2)/total(1.0/toiii_err^2)
;;    out.toiii_err = 1.0/sqrt(total(1.0/toiii_err^2))
;;
;;    out.log12oh = total(log12oh/log12oh_err^2)/total(1.0/log12oh_err^2)
;;    out.log12oh_err = 1.0/sqrt(total(1.0/log12oh_err^2))
;;
