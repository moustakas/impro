pro parse_bresolin02, out
; jm06jan18uofa - parse Table 2 from Bresolin & Kennicutt 2002 into
;                 2002_bresolin.sex 

    outpath = hiiregions_path()

    data = rsex(outpath+'02bresolin/raw_bresolin02.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.hii_galaxy,/remove)
    out.hii_region = strcompress(data.hii_region,/remove)
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset

    good = where(data.oii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727 = data[good].oii/100.0
       out[good].oii_3727_err = 0.05*out[good].oii_3727
    endif
    
    good = where(data.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = data[good].oiii/100.0
       out[good].oiii_5007_err = 0.05*out[good].oiii_5007
    endif
    
    good = where(data.nii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = data[good].nii/100.0
       out[good].nii_6584_err = 0.05*out[good].nii_6584
    endif
    
    good = where(data.sii_6716 gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6716 = data[good].sii_6716/100.0
       out[good].sii_6716_err = 0.05*out[good].sii_6716
    endif
    
    good = where(data.sii_6731 gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6731 = data[good].sii_6731/100.0
       out[good].sii_6731_err = 0.05*out[good].sii_6731
    endif
    
    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '2002_bresolin.sex'
    reference = 'Bresolin & Kennicutt 2002, ApJ, 572, 838'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_BRESOLIN02 '+im_today()+'.'
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
    
