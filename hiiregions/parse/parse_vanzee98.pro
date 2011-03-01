pro parse_vanzee98
; jm04feb10uofa - parse Tables 3 and 5 from van Zee et al. (1998) to
;                 pull out the columns that go into 1998_vanzee.dat 
; jm06jan17uofa - updated
; jm06mar20uofa - fixed some indexing problems

    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
;   oratio = 2.88 ; <-- this is what Lisa uses
    nratio = branch.n_ii  ; 3.054
    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio

    t1 = mrdfits(outpath+'98vanzee/vanZee_table1.fits',1,/silent)
    t3 = mrdfits(outpath+'98vanzee/vanZee_table3.fits',1,/silent)
    t4 = mrdfits(outpath+'98vanzee/vanZee_table4.fits',1,/silent)
    t5 = mrdfits(outpath+'98vanzee/vanZee_table5.fits',1,/silent)
    nobject = n_elements(t3)
    ngalaxy = n_elements(t1)
    
; Tables 3 and 4 are matched but table 5 is missing two objects,
; +053-007 in NGC1637 and -003-003 in NGC4395; also, the coordinates
; of one HII region in NGC1232 are +117+084 in Tables 3 and 4, but
; +117-084 in Table 5; assume Tables 3 and 4 are correct; finally,
; -018-009 and +021+011 are row-switched in Table 5, relative to
; Tables 3 and 4!

    fix1 = where(strmatch(t5._vsh98_,'*+117-084*'),nfix)
    t5[fix1]._vsh98_ = '+117+084'

    fix_t5_1 = where(strmatch(t5._vsh98_,'*+021+011*'))
    fix_t5_2 = where(strmatch(t5._vsh98_,'*-018-009*'))

    tempt5 = t5
    tempt5[fix_t5_1] = t5[fix_t5_2]
    tempt5[fix_t5_2] = t5[fix_t5_1]
    t5 = tempt5
    
;   niceprint, t5.ngc_ic, strtrim(t3._vsh98_,2), strtrim(t5._vsh98_,2), long(strtrim(t5._vsh98_,2) eq strtrim(t3._vsh98_,2))
    missing = cmset_op(strtrim(t3._vsh98_,2),'and',/not2,strtrim(t5._vsh98_,2),/index)
    good = lindgen(nobject)
    remove, missing, good

    newt5 = im_empty_structure(t5,ncopies=nobject)
    newt5[good] = t5
    newt5[missing].recno   = t3[missing].recno
    newt5[missing].ngc_ic  = t3[missing].ngc_ic
    newt5[missing]._vsh98_ = t3[missing]._vsh98_
    newt5[missing].slit    = t3[missing].slit
    newt5[missing].rad     = t4[missing].rad

;   niceprint, newt5.ngc_ic, strtrim(t3._vsh98_,2), strtrim(newt5._vsh98_,2), long(strtrim(newt5._vsh98_,2) eq strtrim(t3._vsh98_,2))
    t5 = newt5
    
; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = strcompress(t3.ngc_ic,/remove)
    out.hii_region = t3._vsh98_

    out.raoffset = strmid(out.hii_region,0,4)
    out.deoffset = strmid(out.hii_region,4)
;   niceprint, out.hii_galaxy, out.hii_region, out.raoffset, out.deoffset

; R/R0

    for i = 0L, ngalaxy-1L do begin

       match = where(t1[i].ngc_ic eq t5.ngc_ic,nmatch)
;      print, t1[i].ngc_ic, nmatch

       good = where(t1[i].r25 gt 0.0,ngood)
       if ngood ne 0L then out[match].radius = t5[match].rad/60.0
       if ngood ne 0L then out[match].rr25 = t5[match].rad/t1[i].r25

    endfor

; [O II] 3727    
    
    g = where(t3._oii_ gt 0.0)
    out[g].oii_3727 = t3[g]._oii_
    out[g].oii_3727_err = t3[g].e__oii_

; [O III] 4363

    g = where(t4._oiii_ gt 0.0)
    out[g].oiii_4363 = t3[g]._oiii_/t4[g]._oiii_
    out[g].oiii_4363_err = im_compute_error(t3[g]._oiii_,t3[g].e__oiii_,t4[g]._oiii_,t4[g].e__oiii_,/quotient)

; [O III] 5007    
    
    g = where(t3._oiii_ gt 0.0)
    oiii = t3[g]._oiii_ ; 4959 + 5007
    oiii_err = t3[g].e__oiii_
    star = where(t3[g].n__oiii_ eq '*',comp=nostar)

    oiii5007 = oiii
    oiii5007_err = oiii_err
    
    oiii5007[star] = oiii[star]/ocor
    oiii5007_err[star] = oiii_err[star]/ocor/sqrt(2.0) ; note factor of 1.4
    oiii5007[nostar] = oiii[nostar]/ocor
    oiii5007_err[nostar] = oiii_err[nostar]/ocor/sqrt(2.0) ; note factor of 1.4

    out[g].oiii_5007 = oiii5007
    out[g].oiii_5007_err = oiii5007_err

; [N II] 6584
    
    g = where(t3._nii_ gt 0.0)
    out[g].nii_6584 = t3[g]._nii_/ncor
    out[g].nii_6584_err = t3[g].e__nii_/ncor/sqrt(2.0) ; note factor of 1.4

; [S II] 6716, 6731

    g = where(t4._sii_ gt 0.0,ng)
    cor6716 = t4[g]._sii_/(1+t4[g]._sii_)
    cor6716_err = im_compute_error(t4[g]._sii_,t4[g].e__sii_,1+t4[g]._sii_,t4[g].e__sii_,/quotient)

    cor6731 = 1.0/(1+t4[g]._sii_)
    cor6731_err = im_compute_error(fltarr(ng)+1.0,fltarr(ng),1+t4[g]._sii_,t4[g].e__sii_,/quotient)

    out[g].sii_6716 = t3[g]._sii_*cor6716
    out[g].sii_6716_err = sqrt(t3[g].e__sii_^2 + cor6716_err^2)

    out[g].sii_6731 = t3[g]._sii_*cor6731
    out[g].sii_6731_err = sqrt(t3[g].e__sii_^2 + cor6731_err^2)
    
; H-alpha
    
    g = where(t3.halpha gt 0.0)
    out[g].ha = t3[g].halpha
    out[g].ha_err = t3[g].e_halpha

; T_e and abundances determined directly from [O III] 4363; +088-119
; in NGC4395 does not have the uncertainty properly propagated; fix
; this 

    g = where((t4._oiii_ gt 0.0) and (t5.t_o___ gt 0.0))
    niceprint, t4[g].ngc_ic, t4[g]._vsh98_, t5[g]._vsh98_, t5[g].t_o___, t5[g].e_t_o___
    
    fix = where(strmatch(t5[g]._vsh98_,'*+088-119*'),nfix)
    t5[g[fix]].e_t_o___ = 3495.0

    out[g].t4363 = t5[g].t_o___/1E4
    out[g].t4363_err = t5[g].e_t_o___/1E4

;   noerr = where(float(t5[g].e_t_o___) eq 0.0)
;   out[g[noerr]].toiii_err = 0.10*out[g[noerr]].toiii ; 10%
    
    out[g].log12oh = t5[g]._12_log_o_h_
    out[g].log12oh_err = t5[g].e_12_log_o_h_

; EW(Hb)    
    
    g = where(t3.ewhbeta gt 0.0)
    out[g].ewhb = t3[g].ewhbeta*10 ; [Angstrom]
    
; remove one typo and NGC4395 -003-003, which is an AGN nucleus (cf Ho
; et al. 1995)

    keep = lindgen(nobject)
    rem = where((strmatch(out.hii_region,'*-003-003*',/fold)) or (strmatch(out.hii_region,'*+203-041*',/fold)))
    remove, rem, keep
    out = out[keep]
    nobject = n_elements(out)
    
; write out

    filename = '1998_vanzee.sex'
    reference = 'van Zee et al. 1998, AJ, 116, 2805'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_VANZEE98 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
    
