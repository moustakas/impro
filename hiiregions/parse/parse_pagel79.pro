pro parse_pagel79, out
; jm06dec13nyu

    outpath = hiiregions_path()

    data = rsex(outpath+'79pagel/raw_pagel79.sex')
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.hii_galaxy,/remove)
    out.hii_region = strcompress(data.hii_region,/remove)
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset
    out.radius = data.radius
    out.rr25 = data.rr25

    good = where(data.oii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727 = 10^data[good].oii
       out[good].oii_3727_err = alog(10.0)*out[good].oii_3727*data[good].oii_err
    endif

    good = where(data.oiii_4363 gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_4363 = 10^data[good].oiii_4363
       out[good].oiii_4363_err = alog(10.0)*out[good].oiii_4363*data[good].oiii_4363_err
    endif

    good = where(data.ha gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].ha = 10^data[good].ha
       out[good].ha_err = alog(10.0)*out[good].ha*data[good].ha_err
    endif

    good = where(data.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = 10^data[good].oiii
       out[good].oiii_5007_err = alog(10.0)*out[good].oiii_5007*data[good].oiii_err
    endif

    good = where(data.nii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = 10^data[good].nii
       out[good].nii_6584_err = alog(10.0)*out[good].nii_6584*data[good].nii_err
    endif

    good = where(data.sii_6716 gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6716 = 10^data[good].sii_6716
       out[good].sii_6716_err = alog(10.0)*out[good].sii_6716*data[good].sii_6716_err
    endif

    good = where(data.sii_6731 gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6731 = 10^data[good].sii_6731
       out[good].sii_6731_err = alog(10.0)*out[good].sii_6731*data[good].sii_6731_err
    endif

    out.ewhb = data.ewhb

; write out

    filename = '1979_pagel.sex'
    reference = 'Pagel et al. 1979, MNRAS, 189, 95'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_PAGEL79 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    printf, lun, '## Data on NGC1365 not tabulated.'
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
    
