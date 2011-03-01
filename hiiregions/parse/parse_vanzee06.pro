pro parse_vanzee06
; jm06dec14nyu

; do not tabulate the [S II] flux ratios    
    
    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    nratio = branch.n_ii  ; 3.054
    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio

    data = rsex(outpath+'06vanzee/raw_vanzee06.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = data.galaxy
    out.hii_region = data.hii_region
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset
    out.ewhb = data.ewhb

    out.oii_3727 = data.oii
    out.oii_3727_err = data.oii_err

    g = where(data.oiii gt -900.0,ng)
    if (ng ne 0L) then begin
       out[g].oiii_5007 = data[g].oiii/ocor
       out[g].oiii_5007_err = data[g].oiii_err/ocor/sqrt(2.0) ; note factor of 1.4
    endif

    g = where(data.nii gt -900.0,ng)
    if (ng ne 0L) then begin
       out[g].nii_6584 = data[g].nii/ncor
       out[g].nii_6584_err = data[g].nii_err/ncor/sqrt(2.0) ; note factor of 1.4
    endif

    out.ha = data.ha
    out.ha_err = data.ha_err

; [O III] 4363; set the electron temperatures and "direct" abundances
; for objects without 4363 detections to -999.0

    out.t4363 = data.toiii/1E4
    out.t4363_err = data.toiii_err/1E4
    out.log12oh = data.log12oh
    out.log12oh_err = data.log12oh_err
    
    g = where(data.oiii_ratio gt -900.0,ng,comp=bad,ncomp=nbad)
    if (ng ne 0L) then begin
       out[g].oiii_4363 = data[g].oiii/data[g].oiii_ratio
       out[g].oiii_4363_err = im_compute_error(data[g].oiii,data[g].oiii_err,data[g].oiii_ratio,data[g].oiii_ratio_err,/quotient)
    endif

    if (nbad ne 0L) then begin
       out[bad].t4363 = -999.0
       out[bad].t4363_err = -999.0
       out[bad].log12oh = -999.0
       out[bad].log12oh_err = -999.0
    endif

; write out

    filename = '2006_vanzee.sex'
    reference = 'van Zee & Haynes 2006, ApJ, 636, 214'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_VANZEE06 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
    
