pro parse_bresolin99, out
; jm05may09uofa
; parse Table 3 from Bresolin, Kennicutt, & Garnett 1999 into
; 1999_bresolin.dat 

    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    nratio = branch.n_ii  ; 3.054
    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio

    data = mrdfits(outpath+'99bresolin/raw_bresolin99.fits',1,/silent)
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    m31 = where(strmatch(data.name,'*messier031*',/fold) eq 1B) & data[m31].name = 'NGC0224'
    m81 = where(strmatch(data.name,'*messier081*',/fold) eq 1B) & data[m81].name = 'NGC3031'
    m51 = where(strmatch(data.name,'*messier051*',/fold) eq 1B) & data[m51].name = 'NGC5194'
    m33 = where(strmatch(data.name,'*messier033*',/fold) eq 1B) & data[m33].name = 'NGC0598'
    out.hii_galaxy = strcompress(data.name,/remove)

    out.hii_region = strcompress(data.hiiname,/remove)
    out.raoffset = data.raoff
    out.deoffset = data.deoff

    noname = where(strcompress(data.hiiname,/remove) eq '')
    out[noname].hii_region = repstr(string(data[noname].raoff,format='(I4.3)')+string(data[noname].deoff,format='(I4.3)'),' ','+')

; several HII regions in M51 has no [O II] detection, causing a
; divide-by-zero in WRITE_HII_REGIONS    
    
    good = where((data.oii gt 0.0) and (data.e_oii gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].oii_3727 = data[good].oii/100.0
       out[good].oii_3727_err = data[good].e_oii/100.0
    endif
    
; ditto
    
    good = where((data.oiii gt 0.0) and (data.e_oiii gt 0.0),ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = data[good].oiii/ocor/100.0
       out[good].oiii_5007_err = data[good].e_oiii/ocor/100.0/sqrt(2.0) ; note factor of 1.4
    endif
    
    good = where((data.nii gt 0.0) and (data.e_nii gt 0.0),ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = data[good].nii/ncor/100.0
       out[good].nii_6584_err = data[good].e_nii/ncor/100.0/sqrt(2.0) ; note factor of 1.4
    endif
    
    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '1999_bresolin.sex'
    reference = 'Bresolin, Kennicutt, & Garnett 1999, ApJ, 510, 104'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_BRESOLIN99 '+im_today()+'.'
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
    
