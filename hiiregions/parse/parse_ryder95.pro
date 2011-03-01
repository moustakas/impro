pro parse_ryder95
; jm05may10uofa - parse Table 3 from Ryder 1995 to generate
;                 1995_ryder.dat 

; no HII-region coordinates are available!
    
    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
    nratio = branch.n_ii  ; 3.054
    ocor = 1.0+1.0/oratio
    ncor = 1.0+1.0/nratio

    data = rsex(outpath+'95ryder/raw_ryder95.sex')
    nobject = n_elements(data)

;   readcol, 'raw_ryder95.dat', galaxy, hiiregion, oiihb, oiihb_err, $
;     oiiihb, oiiihb_err, niihb, niihb_err, hahb, hahb_err, rr25, r25, $
;     /silent, format='A,A,F,F,F,F,F,F,F,F,F,F', comment='#'
;   nobject = n_elements(hiiregion)

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = data.hii_galaxy
    out.hii_region = data.hii_region
    out.radius = data.radius/60.0 ; [arcmin]
    out.rr25 = data.rr25

; the factor of 1.27 converts pixels to arcseconds; first rotate the
; whole coordinate system by -90 degrees; the definition of the
; coordinate axes is a little funny here, but is correct as is,
; including the apparently switched raoffset and deoffset (see
; RYDERTEST) (jm06mar31uofa)

    xy = im_offset_and_rotate(transpose([[data.x-data.xnuc],[data.y-data.ynuc]]),-90.0)
    
    out.deoffset = reform(xy[0,*])
    out.raoffset = reform(xy[1,*])

;   indx = lindgen(10)
;   niceprint, out[indx].hii_region, out[indx].raoffset, out[indx].deoffset

; verify the deprojection code

;   w = where(strmatch(out.hii_galaxy,'*2442*'),nhii)
;   incl = 62.0
;   pa = 90.0+112.0 ; <-- note!
;
;   radius = im_hiiregion_deproject(incl,pa,out[w].raoffset,out[w].deoffset,hii_phi=phi)
;   niceprint, out[w].hii_galaxy, out[w].hii_region, data[w].radius, radius, data[w].hii_phi, phi

    good = where(data.oii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oii_3727 = data[good].oii
       out[good].oii_3727_err = data[good].oii_err
    endif

    good = where(data.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].oiii_5007 = data[good].oiii/ocor
       out[good].oiii_5007_err = data[good].oiii_err/ocor/sqrt(2.0) ; note factor of 1.4
    endif

    good = where(data.nii gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].nii_6584 = data[good].nii/ncor
       out[good].nii_6584_err = data[good].nii_err/ncor/sqrt(2.0) ; note factor of 1.4
    endif

    good = where(data.hahb gt -900,ngood)
    if (ngood ne 0L) then begin
       out[good].ha = data[good].hahb
       out[good].ha_err = data[good].hahb_err
    endif

; write out

    filename = '1995_ryder.sex'
    reference = 'Ryder 1995, ApJ, 444, 610'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_RYDER95 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments

return
end
    
