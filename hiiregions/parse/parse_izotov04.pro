pro parse_izotov04
; jm05may12uofa
; parse Tables 3 and 4 from Izotov et al. 2004 into 2004_izotov.dat 

    outpath = hiiregions_path()

    nobject = 33L
    out = init_hii_region_structure(nobject)

    wave = [3727,4363,5007,6563,6583,6717,6731]
    suffix = ['a','b','c','d','e','f','g','h','i']

    readcol, outpath+'04izotov/raw_izotov04_extradata.dat', galaxy, region, $
      te, te_err, oh, oh_err, format='A,A,F,F,F,F', comment='#', /silent

    indx = 0L
    for i = 0L, n_elements(suffix)-1L do begin

       table1 = im_read_fmr(outpath+'04izotov/raw_izotov04_table3'+suffix[i]+'.dat')
       keep = cmset_op(table1.wave,'AND',wave,/index)
       tags = where(strmatch(tag_names(table1[0]),'F*') eq 1B,nobj)
       ewtags = where(strmatch(tag_names(table1[0]),'EW*') eq 1B)

       hbindx = where(table1.wave eq 4861)
       
       for j = 0L, nobj-1L do begin

          out[indx].hii_galaxy = galaxy[indx]
          out[indx].hii_region = region[indx]
          out[indx].t4363 = te[indx]/1E4
          out[indx].t4363_err = te_err[indx]/1E4
          out[indx].log12oh = oh[indx]
          out[indx].log12oh_err = oh_err[indx]

; line fluxes          
          
          if (table1[keep[0]].(tags[j]) gt -900.0) then begin
             out[indx].oii_3727     = table1[keep[0]].(tags[j])/100.0
             out[indx].oii_3727_err = table1[keep[0]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[1]].(tags[j]) gt -900.0) then begin
             out[indx].oiii_4363     = table1[keep[1]].(tags[j])/100.0
             out[indx].oiii_4363_err = table1[keep[1]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[2]].(tags[j]) gt -900.0) then begin
             out[indx].oiii_5007     = table1[keep[2]].(tags[j])/100.0
             out[indx].oiii_5007_err = table1[keep[2]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[3]].(tags[j]) gt -900.0) then begin
             out[indx].ha     = table1[keep[3]].(tags[j])/100.0
             out[indx].ha_err = table1[keep[3]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[4]].(tags[j]) gt -900.0) then begin
             out[indx].nii_6584     = table1[keep[4]].(tags[j])/100.0
             out[indx].nii_6584_err = table1[keep[4]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[5]].(tags[j]) gt -900.0) then begin
             out[indx].sii_6716     = table1[keep[5]].(tags[j])/100.0
             out[indx].sii_6716_err = table1[keep[5]].(tags[j]+1L)/100.0
          endif
          if (table1[keep[6]].(tags[j]) gt -900.0) then begin
             out[indx].sii_6731     = table1[keep[6]].(tags[j])/100.0
             out[indx].sii_6731_err = table1[keep[6]].(tags[j]+1L)/100.0
          endif

; EW(Hb)

          if (table1[hbindx].(ewtags[j]) gt -900) then out[indx].ewhb = table1[hbindx].(ewtags[j])

          indx = indx + 1L

       endfor
       
    endfor

; write out

    filename = '2004_izotov.sex'
    reference = 'Izotov & Thuan 2004, ApJ, 500, 188'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_IZOTOV04 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments
    
return
end
