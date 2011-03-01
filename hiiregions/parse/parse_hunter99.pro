pro parse_hunter99, out
; jm06jan18uofa - parse Table 3 from Hunter & Hoffman 1999

    splog, 'DO NOT USE -- not many [OII] measurements and fairly noisy data.'
    
    outpath = hiiregions_path()

    data = rsex(outpath+'99hunter/raw_hunter99.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.hii_galaxy,/remove)
    out.hii_region = strcompress(data.hii_region,/remove)

    out.raoffset = data.raoffset
    out.deoffset = data.deoffset
    
    good = where(data.oii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727 = data[good].oii
       out[good].oii_3727_err = data[good].oii_err
    endif
    
    good = where(data.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = data[good].oiii
       out[good].oiii_5007_err = data[good].oiii_err
    endif
    
    good = where(data.nii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = data[good].nii
       out[good].nii_6584_err = data[good].nii_err
    endif
    
    out.ewhb = data.ewhb

    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '1999_hunter.sex'
    reference = 'Hunter & Hoffman 1999 AJ, 117, 2789'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_HUNTER99 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   comments = '# '
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments
    
return
end
    
