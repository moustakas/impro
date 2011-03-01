pro parse_izotov06
; jm07nov02nyu - parse Izotov et al. 2006

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984

    outpath = hiiregions_path()
    table = im_read_vizier_tsv(outpath+'06izotov/2006_izotov.dat')

    nobject = n_elements(table)
    out = init_hii_region_structure(nobject)

; read the NED galaxy names

    nedalldata = djs_readlines(outpath+'06izotov/2006_izotov.ned')
    nned = n_elements(nedalldata)
    ned = replicate({galaxy: '', ra: 0.0D, dec: 0.0D},nned)
    for ii = 0L, nned-1L do begin
       nedall = strsplit(nedalldata[ii],'|',/extract)
       ned[ii].galaxy = strcompress(nedall[1],/remove)
       ned[ii].ra = nedall[2]
       ned[ii].dec = nedall[3]
    endfor

; now cross-match the coordinates to figure out which NED names belong

    splog, 'Searching the distance catalog.'
    ntot = im_djs_angle_match(table._raj2000,table._dej2000,ned.ra,ned.dec,$
      dtheta=10.0/3600.0,units='degrees',mcount=mcount,$
      mindx=mindx,mdist=mdist,mmax=1)
    splog, 'Matched ', ntot, ' galaxies against the NED catalog.'
;   if (ntot ne ngalaxy) then message, 'Problem here.'
    
    good = where(mindx ne -1L,ngood,comp=bad)
    srt = sort(mdist[good])
;   niceprint, table[good[srt]].sdss, ned[mindx[good[srt]]].galaxy, mdist[good[srt]]*3600.0

    out[good].hii_galaxy  = ned[mindx[good]].galaxy
    out[good].hii_region  = 'SDSS'+strtrim(table[good].sdss,2)
    out[good].t4363       = table[good].te
    out[good].t4363_err   = table[good].e_te
    out[good].log12oh     = table[good]._12_logo_h
    out[good].log12oh_err = table[good].e_12_logo_h

; not everything matches, however; the following objects are
; individual HII regions in M101=NGC5457:
;   SDSSJ140301.17+541429.3=NGC5457:[HK83]315=NGC5455
;   SDSSJ140332.35+541719.4=NGC5457:[HK83]060
;   SDSSJ140429.48+542347.2=NGC5471 

;   niceprint, table[bad].sdss, table[bad].name
    out[bad].hii_galaxy  = 'MESSIER101'
    out[bad].hii_region  =  'SDSS'+strtrim(table[bad].sdss,2)
    out[bad].t4363       = table[bad].te
    out[bad].t4363_err   = table[bad].e_te
    out[bad].log12oh     = table[bad]._12_logo_h
    out[bad].log12oh_err = table[bad].e_12_logo_h
    
; line-fluxes

    good = where(table.f3727 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727     = table[good].f3727
       out[good].oii_3727_err = table[good].e_f3727
    endif
    
    good = where(table.f7325 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_7325     = table[good].f7325
       out[good].oii_7325_err = table[good].e_f7325
    endif
    
    good = where(table.f4363 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_4363     = table[good].f4363
       out[good].oiii_4363_err = table[good].e_f4363
    endif
    
    good = where(table.f4959 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007     = table[good].f4959*(1.0+oratio)
       out[good].oiii_5007_err = table[good].e_f4959*(1.0+oratio)
    endif
    
    good = where(table.f6563 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].ha     = table[good].f6563
       out[good].ha_err = table[good].e_f6563
    endif
    
    good = where(table.f6584 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584     = table[good].f6584
       out[good].nii_6584_err = table[good].e_f6584
    endif

    good = where(table.f6312 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].siii_6312     = table[good].f6312
       out[good].siii_6312_err = table[good].e_f6312
    endif

    good = where(table.f9069 gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].siii_9069     = table[good].f9069
       out[good].siii_9069_err = table[good].e_f9069
    endif

; [S II] lines not measured

;   good = where(table.f6716 gt 0.0,ngood)
;   if (ngood ne 0L) then begin
;      out[good].sii_6716     = table[good].f6716
;      out[good].sii_6716_err = table[good].e_f6716
;   endif
;   
;   good = where(table.f6731 gt 0.0,ngood)
;   if (ngood ne 0L) then begin
;      out[good].sii_6731     = table[good].f6731
;      out[good].sii_6731_err = table[good].e_f6731
;   endif

    good = where(table.ewhb gt 0.0,ngood)
    if (ngood ne 0L) then begin
       out[good].ewhb = table[good].ewhb
    endif

; write out

    filename = '2006_izotov.sex'
    reference = 'Izotov et al. 2006, A&A, 448, 955'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_IZOTOV06 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
