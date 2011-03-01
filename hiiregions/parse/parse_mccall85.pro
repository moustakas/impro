pro parse_mccall85
; jm04jun25uofa - parse Table 2 from McCall et al. 1985 into
;                 1985_mccall.dat 
; jm06jan17uofa - updated

    outpath = hiiregions_path()

    data = rsex(outpath+'85mccall/raw_mccall85.sex')
    nobject = n_elements(data)

    log_oii = data.oii
    log_oii_err = data.oii_err
    log_oiii4363 = data.oiii_4363
    log_oiii4363_err = data.oiii_4363_err
    log_oiii5007 = data.oiii
    log_oiii5007_err = data.oiii_err
    log_ha = data.ha
    log_ha_err = data.ha_err
    log_nii = data.nii
    log_nii_err = data.nii_err
    log_sii6716 = data.sii_6716
    log_sii6716_err = data.sii_6716_err
    log_sii6731 = data.sii_6731
    log_sii6731_err = data.sii_6731_err

;   readcol, 'raw_mccall85.dat', galaxy, hiiregion, log_oii, log_oii_err, $
;     log_oiii4363, log_oiii4363_err, log_oiii5007, log_oiii5007_err, $
;     log_ha, log_ha_err, log_nii, log_nii_err, log_sii6716, log_sii6716_err, $
;     log_sii6731, log_sii6731_err, rr0, ewhb, radius, $
;     format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', $
;     comment='#', /silent
;   nobject = n_elements(galaxy)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = data.galaxy
    out.hii_region = repstr(string(data.raoffset,format='(I4.3)')+string(data.deoffset,format='(I4.3)'),' ','+')
    star = where(strmatch(out.hii_region,'-606*'),nstar)
    if (nstar ne 0L) then out[star].hii_region = repstr(string(data[star].raoffset,format='(I5.4)')+$
      string(data[star].deoffset,format='(I5.4)'),' ','+')
    
    out.ewhb = data.ewhb
    out.radius = data.radius ; [arcmin]
    out.rr25 = data.rr25
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset

    good = where(log_oii gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727 = 10^log_oii[good]/100.0
       out[good].oii_3727_err = alog(10.0)*log_oii_err[good]*out[good].oii_3727
    endif

    good = where(log_oiii4363 gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_4363 = 10^log_oiii4363[good]/100.0
       out[good].oiii_4363_err = alog(10.0)*log_oiii4363_err[good]*out[good].oiii_4363
    endif

    good = where(log_oiii5007 gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = 10^log_oiii5007[good]/100.0
       out[good].oiii_5007_err = alog(10.0)*log_oiii5007_err[good]*out[good].oiii_5007
    endif

    good = where(log_nii gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = 10^log_nii[good]/100.0
       out[good].nii_6584_err = alog(10.0)*log_nii_err[good]*out[good].nii_6584
    endif

    good = where(log_sii6716 gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6716 = 10^log_sii6716[good]/100.0
       out[good].sii_6716_err = alog(10.0)*log_sii6716_err[good]*out[good].sii_6716
    endif

    good = where(log_sii6731 gt -900.0,ngood)
    if (ngood ne 0L) then begin
       out[good].sii_6731 = 10^log_sii6731[good]/100.0
       out[good].sii_6731_err = alog(10.0)*log_sii6731_err[good]*out[good].sii_6731
    endif

    good = where(log_ha gt -900.0,ngood,comp=bad,ncomp=nbad)
    if (ngood ne 0L) then begin
       out[good].ha = 10^log_ha[good]/100.0
       out[good].ha_err = alog(10.0)*log_ha_err[good]*out[good].ha
    endif
    if (nbad ne 0L) then begin
       out[bad].ha = 2.86 ; (return_tbalmer(/HaHb))[0]
       out[bad].ha_err = 0.05*out[bad].ha
    endif

;   out.ha = return_tbalmer(/HaHb)
;   out.ha_err = 0.05*out.ha
    
; remove crummy HII regions
;
;   keep = lindgen(nobject)
;   rem = where((strmatch(out.hii_galaxy,'*NGC4395*',/fold) eq 1B) and $
;     (strmatch(out.hii_region,'*H1*',/fold) eq 1B))
;   remove, rem, keep
;   out = out[keep]
;   nobject = n_elements(out)
    
; write out

    filename = '1985_mccall.sex'
    reference = 'Mccall, Rybski, & Shields 1985, ApJS, 57, 1'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_MCCALL85 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments
    
return
end
    
