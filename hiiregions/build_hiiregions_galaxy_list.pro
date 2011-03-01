pro build_hiiregions_galaxy_list, clobber=clobber
; jm10feb22ucsd - build the unique list of galaxies in the current
;   version of the database

    hiipath = hiiregions_path()
    version = hiiregions_version()

    outfile = hiipath+'hiiregions_galaxy_list_'+version+'.txt'
    if file_test(outfile,/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

; read the     
    readcol, hiipath+'datafiles_hii_regions.txt', hii_datafile, $
      format='A', comment='#', delimiter='|', /silent
    readcol, hiipath+'datafiles_hii_galaxies.txt', galaxy_datafile, $
      format='A', comment='#', delimiter='|', /silent

    datafiles = strtrim([hii_datafile,galaxy_datafile],2)
    nfile = n_elements(datafiles)

    for ii = 0, nfile-1 do begin
       info1 = rsex(hiipath+datafiles[ii])
       info1 = struct_trimtags(info1,select='hii_galaxy')
       if (ii eq 0) then info = info1 else info = [info,info1]
    endfor

    gal = strupcase(strtrim(info.hii_galaxy,2))
    ugal = gal[uniq(gal,sort(gal))]

    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '# HII-region database'
    printf, lun, '# Version '+version
    printf, lun, '# Unique list of galaxy names built by BUILD_HIIREGIONS_GALAXY_LIST'
    printf, lun, '# '
    niceprintf, lun, ugal
    free_lun, lun

return
end
    
