pro parse_guseva07
; jm07nov02nyu - parse Tables 2 and 3 from Guseva et al. 2007

    outpath = hiiregions_path()

    table2 = im_read_vizier_tsv(outpath+'07guseva/2007_guseva_table2.dat')
    table3 = im_read_vizier_tsv(outpath+'07guseva/2007_guseva_table3b.dat')

    ngalaxy = n_elements(table2)
    out = init_hii_region_structure(ngalaxy)

; parse the line-flux structure

    nline = n_elements(table3)/ngalaxy
    
; read the NED galaxy names

    nedalldata = djs_readlines(outpath+'07guseva/2007_guseva.ned')
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
    ntot = im_djs_angle_match(table2._raj2000,table2._dej2000,ned.ra,ned.dec,$
      dtheta=75.0/3600.0,units='degrees',mcount=mcount,$
      mindx=mindx,mdist=mdist,mmax=1)
    splog, 'Matched ', ntot, ' galaxies against the NED catalog.'
    if (ntot ne ngalaxy) then message, 'Problem here.'
    
    good = where(mindx ne -1L,ngood)
    srt = sort(mdist[good])
;   niceprint, table2[good[srt]].name, ned[mindx[good[srt]]].galaxy, mdist[good[srt]]*3600.0

    out[good].hii_galaxy  = ned[mindx[good]].galaxy
    out[good].hii_region  = table2[good].name
    out[good].t4363       = table2[good].toiii/1D4
    out[good].t4363_err   = table2[good].e_toiii/1D4
    out[good].log12oh     = table2[good]._12o_h
    out[good].log12oh_err = table2[good].e_12o_h
    
; line-fluxes

    for ig = 0L, ngalaxy-1L do begin
       line = table3[ig*nline:(ig+1L)*nline-1L]
       match = where(line.lambda eq 3727)
       if (line[match].fluxc gt 0.0) and (line[match].e_fluxc gt 0.0) then begin
          out[ig].oii_3727     = line[match].fluxc/100.0
          out[ig].oii_3727_err = line[match].e_fluxc/100.0
       endif
       match = where(line.lambda eq 4363)
       if (line[match].fluxc gt 0.0) and (line[match].e_fluxc gt 0.0) then begin
          out[ig].oiii_4363     = line[match].fluxc/100.0
          out[ig].oiii_4363_err = line[match].e_fluxc/100.0
       endif
       match = where(line.lambda eq 5007)
       if (line[match].fluxc gt 0.0) and (line[match].e_fluxc gt 0.0) then begin
          out[ig].oiii_5007     = line[match].fluxc/100.0
          out[ig].oiii_5007_err = line[match].e_fluxc/100.0
       endif
; assign the H-alpha flux based on the derive Te value and the error
; based on the H-beta flux error
       match = where(line.lambda eq 4861)
       if (line[match].fluxc gt 0.0) and (line[match].e_fluxc gt 0.0) and (line[match].ew gt 0.0) then begin
          out[ig].ha     = return_tbalmer(out[ig].t4363*1D4)
          out[ig].ha_err = out[ig].ha/(line[match].fluxc/line[match].e_fluxc)
          out[ig].ewhb   = line[match].ew
       endif
    endfor

; write out

    filename = '2007_guseva.sex'
    reference = 'Guseva et al. 2007, A&A, 464, 885'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_GUSEVA07 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
