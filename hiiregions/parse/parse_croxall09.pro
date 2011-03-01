pro parse_croxall09, out
; jm10mar08ucsd - parse Croxall et al. 2009; Garland has been removed,
; as have all the supernova remnants

    outpath = hiiregions_path()

    table2 = rsex(outpath+'09croxall/table2.sex')
    table3 = rsex(outpath+'09croxall/table3.sex')
    table4 = rsex(outpath+'09croxall/table4.sex')
    table5 = rsex(outpath+'09croxall/table5.sex')
    nobj = n_elements(table2)

; initialize the output data structure
    out = init_one_hiiregion_structure(nobj)

    out.hii_galaxy = strcompress(table2.galaxy,/remove)
    out.hii_region = strcompress(table2.region,/remove)

; compute the RA/DEC offsets
    gal = [$
      'UGC04305',$ ; =HoII
      'UGC04459',$ ; =DDO53
      'UGC05139',$ ; =HoI
      'KDG061',$   
      'UGC05336',$ ; =HoIX
      'UGC05666',$ ; =IC2574
      'UGC05692',$ ; =DDO82
      'UGC05918',$ ; =DDO87
      'UGC08201']  ; =DDO165
    ra = [$
      '08:19:05.0',$
      '08:34:07.2',$
      '09:40:32.3',$
      '09:57:03.1',$
      '09:57:32.0',$
      '10:28:23.5',$
      '10:30:35.0',$
      '10:49:36.5',$
      '13:06:24.8']
    dec = [$
      '+70:43:12',$
      '+66:10:54',$
      '+71:10:56',$
      '+68:35:31',$
      '+69:02:45',$
      '+68:24:44',$
      '+70:37:07',$
      '+65:31:50',$
      '+67:42:25']
    
    for ii = 0, nobj-1 do begin
       ww = where(strtrim(out[ii].hii_galaxy,2) eq gal)
       cosd = cos(im_hms2dec(dec[ww])*!dtor)
       out[ii].deoffset = (im_hms2dec(table2[ii].dec)-im_hms2dec(dec[ww]))*3600.0
       out[ii].raoffset = (im_hms2dec(table2[ii].ra)-im_hms2dec(ra[ww]))*15.0*cosd*3600.0
    endfor

; fill the table
    struct_assign, table3, out, /nozero
    struct_assign, table4, out, /nozero
    struct_assign, table5, out, /nozero

; T(4363) needs special care; if T(4363) was not detected then
; don't use the electron temperature or oxygen abundance from Table 5 
    out.t4363 = out.t4363/1E4
    for ii = 0, nobj-1 do out[ii].t4363_err = $
      mean([table5[ii].t4363_up_err,table5[ii].t4363_lo_err])/1E4

    toss = where(out.oiii_4363 lt -900)
    out[toss].t4363 = -999.0
    out[toss].t4363_err = -999.0
    out[toss].log12oh = -999.0
    out[toss].log12oh_err = -999.0
    
; write out
    filename = '2009_croxall.sex'
    reference = 'Croxall et al. 2009, ApJ, 705, 723'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_CROXALL09 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
    
