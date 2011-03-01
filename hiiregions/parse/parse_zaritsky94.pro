pro parse_zaritsky94
; jm04feb9uofa - parse Table 2 from Zaritsky, Kennicutt, & Huchra 1994
;                into 1994_zaritsky.dat 
; jm06jan17uofa - updated

    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    ocor = 1.0+1.0/oratio

    data = rsex(outpath+'94zaritsky/raw_zaritsky94.sex')
    nobject = n_elements(data)

;   readcol, 'raw_zaritsky94.dat', galaxy, oiihb, oiihb_err, $
;     oiiihb, oiiihb_err, rr0, sii6716, sii6716_err, sii6731, $
;     sii6731_err, r0, region, format='A,D,D,D,D,D,D,D,D,D,F,A', $
;     comment='#', /silent
;   nobject = n_elements(galaxy)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = data.galaxy
    out.hii_region = repstr(string(data.raoffset,format='(I4.3)')+string(data.deoffset,format='(I4.3)'),' ','+')
    out.radius = data.rr25*data.r25
    out.rr25 = data.rr25
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset

    out.oii_3727 = data.oii
    out.oii_3727_err = data.oii_err

    out.oiii_5007 = data.oiii/ocor
    out.oiii_5007_err = data.oiii_err/ocor/sqrt(2.0) ; note factor of 1.4

    out.sii_6716 = data.sii_6716
    out.sii_6716_err = data.sii_6716_err

    out.sii_6731 = data.sii_6731
    out.sii_6731_err = data.sii_6731_err

    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '1994_zaritsky.sex'
    reference = 'Zaritsky, Kennicutt, & Huchra 1994, ApJ, 420, 87'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_ZARITSKY94 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments

return
end
    
