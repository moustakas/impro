pro parse_kennicutt96
; jm04jun24uofa - parse Table 2 from Kennicutt & Garnett 1996 to
;                 generate 1996_kennicutt.dat 

    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    nratio = branch.n_ii  ; 3.054
    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio

    data = rsex(outpath+'96kennicutt/raw_kennicutt96.sex')
    nobject = n_elements(data)

;   readcol, 'raw_kennicutt96.dat', hiiregion, oiihb, oiihb_err, $
;     oiiihb, oiiihb_err, rr0, niihb, niihb_err, ewhb, r0, /silent, $
;     format='A,F,F,F,F,F,F,F,F,F', comment='#'
;   nobject = n_elements(hiiregion)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = 'NGC5457'
    out.hii_region = data.hii_region
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset
    out.radius = data.radius
    out.rr25 = data.rr25
    out.ewhb = data.ewhb

    out.oii_3727 = data.oii
    out.oii_3727_err = data.oii_err

    out.oiii_5007 = data.oiii/ocor
    out.oiii_5007_err = data.oiii_err/ocor/sqrt(2.0) ; note factor of 1.4

    out.nii_6584 = data.nii/ncor
    out.nii_6584_err = data.nii_err/ncor/sqrt(2.0) ; note factor of 1.4

    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '1996_kennicutt.sex'
    reference = 'Kennicutt & Garnett 1996, ApJ, 456, 504'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_KENNICUTT96 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments

return
end
    
