pro parse_oey93
; jm04jun24uofa - parse Table 2 from Oey & Kennicutt 1993 into
;                 1993_oey.dat  
; jm06jan17uofa - updated

    outpath = hiiregions_path()

    branch = im_branch_ratios()
    oratio = branch.o_iii ; 2.984
;   oratio = 2.98507 ; this is what Sally uses
    ocor = 1.0+1.0/oratio

    data = rsex(outpath+'93oey/raw_oey93.sex')
    nobject = n_elements(data)

    log_oiihb = data.oii
    log_oiihb_err = data.oii_err
    log_oiiihb = data.oiii
    log_oiiihb_err = data.oiii_err
    
;   readcol, 'raw_data93.dat', galaxy, log_oiihb, log_oiihb_err, $
;     log_oiiihb, log_oiiihb_err, rr0, ewhb, logd0, region, $
;     format='A,F,F,F,F,F,F,F,A', comment='#', /silent
;   nobject = n_elements(data)

    oiihb = 10^log_oiihb
    oiihb_err = alog(10.0)*oiihb*log_oiihb_err
    
    oiiihb = 10^log_oiiihb
    oiiihb_err = alog(10.0)*oiiihb*log_oiiihb_err

; initialize the output data structure

    out = init_hii_region_structure(nobject)

    out.hii_galaxy = data.galaxy
    out.hii_region = data.hii_region
    out.raoffset = data.raoffset
    out.deoffset = data.deoffset

    out.rr25 = data.rr25
    out.radius = data.rr25 * (0.1*10^data.logd0/2.0)

    out.ewhb = data.ewhb
    out.oii_3727 = oiihb
    out.oii_3727_err = oiihb_err

    out.oiii_5007 = oiiihb/ocor
    out.oiii_5007_err = oiiihb_err/ocor/sqrt(2.0) ; note factor of 1.4

    out.ha = 2.86 ; return_tbalmer(/hahb)
    out.ha_err = 0.05*out.ha

; write out

    filename = '1993_oey.sex'
    reference = 'Oey & Kennicutt 1993, ApJ, 411, 137'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_OEY93 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
;   write_parsed_hii_region_structure, out, filename, reference, comments=comments
    
return
end
    
