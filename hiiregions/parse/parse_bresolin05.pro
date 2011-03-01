pro parse_bresolin05, out
; jm06jan18uofa - parse Bresolin et al. 2005

    outpath = hiiregions_path()

    n1232 = mrdfits(outpath+'05bresolin/05bresolin_ngc1232.fits.gz',1,/silent)
    n1365 = mrdfits(outpath+'05bresolin/05bresolin_ngc1365.fits.gz',1,/silent)
    n2903 = mrdfits(outpath+'05bresolin/05bresolin_ngc2903.fits.gz',1,/silent)
    n2997 = mrdfits(outpath+'05bresolin/05bresolin_ngc2997.fits.gz',1,/silent)
    n5236 = mrdfits(outpath+'05bresolin/05bresolin_ngc5236.fits.gz',1,/silent)
    sex = rsex(outpath+'05bresolin/raw_bresolin05.sex')
    sex2 = rsex(outpath+'05bresolin/raw_bresolin05_te.sex')
    data = struct_addtags(sex,struct_append(struct_append(struct_append($
      struct_append(n1232,n1365),n2903),n2997),n5236))
    
    nobject = n_elements(data)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(data.hii_galaxy,/remove)
    out.hii_region = strcompress(data.hii_region,/remove)
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset
    out.radius = data.radius/60.0 ; [arcmin]
    out.rr25 = data.radius/data.r25
    out.ewhb = data.ewhb

    for i = 0L, n_elements(sex2)-1L do begin
       match = where((strmatch(out.hii_galaxy,sex2[i].hii_galaxy) eq 1B) and $
         (strmatch(out.hii_region,sex2[i].hii_region) eq 1B),nmatch)
       if (nmatch ne 0L) then begin
; measured, auroral temperatures
          if (sex2[i].t5755 gt -900.0) then begin
             out[match].t5755 = sex2[i].t5755/1E4
             out[match].t5755_err = sex2[i].t5755_err/1E4
          endif
          if (sex2[i].t6312 gt -900.0) then begin
             out[match].t6312 = sex2[i].t6312/1E4
             out[match].t6312_err = sex2[i].t6312_err/1E4
          endif
          if (sex2[i].t7325 gt -900.0) then begin
             out[match].t7325 = sex2[i].t7325/1E4
             out[match].t7325_err = sex2[i].t7325_err/1E4
          endif
; final adopted electron temperatures          
;         out[match].toiii = sex2[i].toiii/1E4
;         out[match].toiii_err = sex2[i].toiii_err/1E4
;         out[match].tsiii = sex2[i].tsiii/1E4
;         out[match].tsiii_err = sex2[i].tsiii_err/1E4
;         out[match].toii = sex2[i].toii/1E4
;         out[match].toii_err = sex2[i].toii_err/1E4
          out[match].log12oh = sex2[i].log12oh
          out[match].log12oh_err = sex2[i].log12oh_err
       endif
    endfor

    good = where((data.f3727 gt 0.0) and (data.e_f3727 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].oii_3727 = data[good].f3727/100.0
       out[good].oii_3727_err = data[good].e_f3727/100.0
    endif

    good = where((data.f7325 gt 0.0) and (data.e_f7325 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].oii_7325 = data[good].f7325/100.0
       out[good].oii_7325_err = data[good].e_f7325/100.0
    endif

    good = where((data.f5007 gt 0.0) and (data.e_f5007 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = data[good].f5007/100.0
       out[good].oiii_5007_err = data[good].e_f5007/100.0
    endif

    good = where((data.f6583 gt 0.0) and (data.e_f6583 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].nii_6584 = data[good].f6583/100.0
       out[good].nii_6584_err = data[good].e_f6583/100.0
    endif

    good = where((data.f6563 gt 0.0) and (data.e_f6563 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].ha = data[good].f6563/100.0
       out[good].ha_err = data[good].e_f6563/100.0
    endif

    good = where((data.f6716 gt 0.0) and (data.e_f6716 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].sii_6716 = data[good].f6716/100.0
       out[good].sii_6716_err = data[good].e_f6716/100.0
    endif

    good = where((data.f6731 gt 0.0) and (data.e_f6731 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].sii_6731 = data[good].f6731/100.0
       out[good].sii_6731_err = data[good].e_f6731/100.0
    endif

    good = where((data.f5755 gt 0.0) and (data.e_f5755 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].nii_5755 = data[good].f5755/100.0
       out[good].nii_5755_err = data[good].e_f5755/100.0
    endif

    good = where((data.f6312 gt 0.0) and (data.e_f6312 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].siii_6312 = data[good].f6312/100.0
       out[good].siii_6312_err = data[good].e_f6312/100.0
    endif

    good = where((data.f9069 gt 0.0) and (data.e_f9069 gt 0.0),ngood) 
    if (ngood ne 0L) then begin
       out[good].siii_9069 = data[good].f9069/100.0
       out[good].siii_9069_err = data[good].e_f9069/100.0
    endif

; write out

    filename = '2005_bresolin.sex'
    reference = 'Bresolin et al. 2005, A&A, 441, 981'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_BRESOLIN05 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
