pro parse_rosolowsky08
; jm08may07nyu - parse Table 1 from Rosolowsky & Simon 2008

    outpath = hiiregions_path()

    table = rsex(outpath+'08rosolowsky/2008_rosolowsky_table1.sex')
    ngalaxy = n_elements(table)
    out = init_hii_region_structure(ngalaxy)

    out.hii_galaxy  = 'NGC0598' ; M33
    out.hii_region  = table.region
    out.t4363       = table.toiii/1D4
    out.t4363_err   = table.toiii_err/1D4
    out.log12oh     = table.log12oh
    out.log12oh_err = table.log12oh_err

    out.oii_3727      = table.oii
    out.oii_3727_err  = table.oii_err
    out.oiii_4363     = table.oiii_4363
    out.oiii_4363_err = table.oiii_4363_err
    out.oiii_5007     = table.oiii_5007
    out.oiii_5007_err = table.oiii_5007_err

    gra = '01h33m50.9s' & gdec = '+30d39m36s'
    cosd = cos(im_hms2dec(gdec)*!dtor)

    out.deoffset = (im_hms2dec(table.dec)-im_hms2dec(gdec))*3600.0
    out.raoffset = (im_hms2dec(table.ra)-im_hms2dec(gra))*15.0*cosd*3600.0

; write out

    filename = '2008_rosolowsky.sex'
    reference = 'Rosolowsky & Simon 2008, ApJ, 675, 1213'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_ROSOLOWSKY08 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
