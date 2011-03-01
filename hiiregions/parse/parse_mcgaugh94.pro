pro parse_mcgaugh94
; jm05may09uofa
; parse Table 3 from McGaugh 1994 into 1994_mcgaugh.dat; the data need
; to be corrected for reddening and Balmer absorption

    outpath = hiiregions_path()
    path = outpath+'94mcgaugh/'
    
    mcgaugh = rsex(path+'raw_mcgaugh94.dat')
    nobject = n_elements(mcgaugh)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = mcgaugh.galaxy
    out.hii_region = mcgaugh.region

; figure out the Balmer absorption correction

    abscor = fltarr(nobject)+1.0
    ew_abs = 2.0

    good = where((mcgaugh.ew_hbeta gt -90.0),ngood)
    if (ngood ne 0L) then begin

       abscor[good] = (1.0 + ew_abs/mcgaugh[good].ew_hbeta)^(-1.0)
       out[good].ewhb = mcgaugh[good].ew_hbeta + ew_abs
    
    endif

; figure out the reddening and compare with McGaugh 

    ebv = fltarr(nobject)
    ebv_err = fltarr(nobject)-1.0

    good = where((mcgaugh.ha_err gt -90.0),ngood)
    if (ngood ne 0L) then begin

       dec = mcgaugh[good].ha*abscor[good]
       dec_err = mcgaugh[good].ha_err*abscor[good]

       ebv1 = get_ebv(dec,decrement_err=dec_err,ebv_err=ebv_err1)
       ebv[good] = ebv1 > 0.0
       ebv_err[good] = ebv_err1

    endif

;   good = where((mcgaugh.c gt -90.0) and (ebv_err gt -1.0))
;   niceprint, ebv[good], 0.78*mcgaugh[good].c, ebv_err[good], 0.78*mcgaugh[good].c_err
;   ploterror, ebv[good], 0.78*mcgaugh[good].c, ebv_err[good], 0.78*mcgaugh[good].c_err, ps=4, $
;     xrange=[-0.1,1.5], yrange=[-0.1,1.5], xsty=3, ysty=3
;   oplot, !x.crange, !y.crange

; --------------------------------------------------
; use McGaugh's reddening values
; --------------------------------------------------
    ebv = fltarr(nobject) & ebv_err = fltarr(nobject)-1.0
    good = where((mcgaugh.c gt -90.0),ngood)
    ebv[good] = 0.78*mcgaugh[good].c
    ebv_err[good] = 0.78*mcgaugh[good].c_err
; --------------------------------------------------
    
; now de-redden the line fluxes; ignore the error in the reddening
    
    good = where((mcgaugh.ha_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(6563.0,/odonnell)))
       out[good].ha = abscor[good]*mcgaugh[good].ha*cor
       out[good].ha_err = abscor[good]*mcgaugh[good].ha_err*cor

    endif

    good = where((mcgaugh.oii_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(3727.0,/odonnell)))
       out[good].oii_3727 = abscor[good]*mcgaugh[good].oii*cor
       out[good].oii_3727_err = abscor[good]*mcgaugh[good].oii_err*cor

    endif

    good = where((mcgaugh.oiii_4363_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(4363.0,/odonnell)))
       out[good].oiii_4363 = abscor[good]*mcgaugh[good].oiii_4363*cor
       out[good].oiii_4363_err = abscor[good]*mcgaugh[good].oiii_4363_err*cor

    endif

    good = where((mcgaugh.oiii_5007_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(5007.0,/odonnell)))
       out[good].oiii_5007 = abscor[good]*mcgaugh[good].oiii_5007*cor
       out[good].oiii_5007_err = abscor[good]*mcgaugh[good].oiii_5007_err*cor

    endif

    good = where((mcgaugh.nii_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(6584.0,/odonnell)))
       out[good].nii_6584 = abscor[good]*mcgaugh[good].nii*cor
       out[good].nii_6584_err = abscor[good]*mcgaugh[good].nii_err*cor

    endif

    good = where((mcgaugh.sii_6716_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(6716.0,/odonnell)))
       out[good].sii_6716 = abscor[good]*mcgaugh[good].sii_6716*cor
       out[good].sii_6716_err = abscor[good]*mcgaugh[good].sii_6716_err*cor

    endif

    good = where((mcgaugh.sii_6731_err gt -90.0) and (ebv_err gt -1.0),ngood)
    if (ngood ne 0L) then begin

       cor = 10^(-0.4*ebv[good]*(k_lambda(4861.0,/odonnell)-k_lambda(6731.0,/odonnell)))
       out[good].sii_6731 = abscor[good]*mcgaugh[good].sii_6731*cor
       out[good].sii_6731_err = abscor[good]*mcgaugh[good].sii_6731_err*cor

    endif

; write out

    filename = '1994_mcgaugh.dat'
    reference = 'McGaugh 1994, ApJ, 426, 135'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_MCGAUGH94 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments

return
end
    
