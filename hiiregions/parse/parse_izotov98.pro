pro parse_izotov98
; jm04jun24uofa
; parse Table 3 from Izotov et al. 1998 into 1998_izotov.dat

    outpath = hiiregions_path()

    readcol, outpath+'98izotov/raw_izotov98.dat', galaxy, altgalaxy, oiihb, oiihb_err, oiii4363hb, oiii4363hb_err, $
      oiii5007hb, oiii5007hb_err, hahb, hahb_err, niihb, niihb_err, sii6716hb, $
      sii6716hb_err, sii6731hb, sii6731hb_err, fhb, ewhb, te, te_err, oh, oh_err, $
      /silent, format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', comment='#'
    nobject = n_elements(galaxy)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = galaxy
    out.hii_region = altgalaxy
    out.ewhb = ewhb
    out.fha = fhb * hahb * 1D16 / 1D14
    out.t4363 = te/1D4
    out.t4363_err = te_err/1D4
    out.log12oh = oh
    out.log12oh_err = oh_err

    out.oii_3727 = oiihb
    out.oii_3727_err = oiihb_err

    out.oiii_4363 = oiii4363hb
    out.oiii_4363_err = oiii4363hb_err

    out.oiii_5007 = oiii5007hb
    out.oiii_5007_err = oiii5007hb_err

    out.ha = hahb
    out.ha_err = hahb_err

    out.nii_6584 = niihb
    out.nii_6584_err = niihb_err

    out.sii_6716 = sii6716hb
    out.sii_6716_err = sii6716hb_err

    out.sii_6731 = sii6731hb
    out.sii_6731_err = sii6731hb_err

; write out

    filename = '1998_izotov.sex'
    reference = 'Izotov & Thuan 1998, ApJ, 500, 188'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_IZOTOV98 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments
    
return
end
