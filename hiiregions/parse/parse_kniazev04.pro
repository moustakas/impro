pro parse_kniazev04
; jm05may12uofa
; parse Tables 3 and 4 from Kniazev et al. 2004 into 2004_kniazev.dat 

    outpath = hiiregions_path()

    table1 = im_read_fmr(outpath+'04kniazev/raw_kniazev04_table1.dat')
    table2 = im_read_fmr(outpath+'04kniazev/raw_kniazev04_table2.dat')
    table3 = im_read_fmr(outpath+'04kniazev/raw_kniazev04_table3.dat')
    table4 = im_read_fmr(outpath+'04kniazev/raw_kniazev04_table4.dat')

    match, table1.id, table2.id, indx1, indx2
    table1 = table1[indx1] & table2 = table2[indx2]
    table3 = table3[indx2] & table4 = table4[indx2]
;   niceprint, table1.id, table2.id, table3.id, table4.id
    
    nobject = n_elements(table1)
    out = init_hii_region_structure(nobject)

    out.hii_galaxy = 'SHOC'+strcompress(table1.id,/remove)
;   out.hii_galaxy = 'SHOC'+string(lindgen(nobject)+1L,format='(I3.3)')
;   out.hii_galaxy = table1.sdss

    altname = where(strcompress(table1.comm,/remove) ne '',comp=none)
    out[none].hii_region = 'H1'       ; table1.comm
    out[altname].hii_region = strtrim(table1[altname].comm,2)

; NED conflicts:    
    
;   SDSSJ020108.28+132660.0 could be "SDSS J020108.28+132700.0"
;   SDSSJ081437.20+490260.0 could be "SDSS J081437.20+490300.0"
;   SDSSJ090760.00+503910.0 could be "SDSS J090800.09+503910.1"
;   SDSSJ094333.84+010659.3 could be "SDSS J094333.84+010659.4"
;   SDSSJ173126.64+591150.2 could be "SDSS J173126.54+591150.2"

;   w = where(strmatch(out.hii_galaxy,'SDSSJ020108.28+132660.0',/fold),nw)
;   if (nw ne 1L) then message, 'Problem!' else out[w].hii_galaxy = 'SDSSJ020108.28+132700.0'
;   w = where(strmatch(out.hii_galaxy,'SDSSJ081437.20+490260.0',/fold),nw)
;   if (nw ne 1L) then message, 'Problem!' else out[w].hii_galaxy = 'SDSSJ081437.20+490300.0'
;   w = where(strmatch(out.hii_galaxy,'SDSSJ090760.00+503910.0',/fold),nw)
;   if (nw ne 1L) then message, 'Problem!' else out[w].hii_galaxy = 'SDSSJ090800.09+503910.1'
;   w = where(strmatch(out.hii_galaxy,'SDSSJ094333.84+010659.3',/fold),nw)
;   if (nw ne 1L) then message, 'Problem!' else out[w].hii_galaxy = 'SDSSJ094333.84+010659.4'
;   w = where(strmatch(out.hii_galaxy,'SDSSJ173126.64+591150.2',/fold),nw)
;   if (nw ne 1L) then message, 'Problem!' else out[w].hii_galaxy = 'SDSSJ173126.54+591150.2'

    good = where(table3.o3727 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727     = table3[good].o3727/100.0
       out[good].oii_3727_err = table3[good].e_o3727/100.0
    endif
    
    good = where(table3.o7320_ gt 0.0,ngood) ; o7320_ = 7320+7330
    if (ngood ne 0L) then begin
       out[good].oii_7325     = table3[good].o7320_/100.0
       out[good].oii_7325_err = table3[good].e_o7320_/100.0
    endif
    
    good = where(table3.o4363 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_4363     = table3[good].o4363/100.0
       out[good].oiii_4363_err = table3[good].e_o4363/100.0
    endif
    
    good = where(table3.o5007 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007     = table3[good].o5007/100.0
       out[good].oiii_5007_err = table3[good].e_o5007/100.0
    endif
    
    good = where(table3.ha gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].ha     = table3[good].ha/100.0
       out[good].ha_err = table3[good].e_ha/100.0
    endif
    
    good = where(table3.s6716 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6716     = table3[good].s6716/100.0
       out[good].sii_6716_err = table3[good].e_s6716/100.0
    endif
    
    good = where(table3.s6731 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6731     = table3[good].s6731/100.0
       out[good].sii_6731_err = table3[good].e_s6731/100.0
    endif
    
    good = where(table2.eqwid gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].ewhb = table2[good].eqwid
    endif
    
    good = where(table4.t_oiii gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].t4363 = table4[good].t_oiii/1E4
       out[good].t4363_err = table4[good].e_t_oiii/1E4
    endif
    
    good = where(table4.o_h gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].log12oh = table4[good].o_h
       out[good].log12oh_err = table4[good].e_o_h
    endif

; write out

    filename = '2004_kniazev.sex'
    reference = 'Kniazev et al. 2004, ApJS, 153, 429'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_KNIAZEV04 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
