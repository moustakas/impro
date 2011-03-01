pro compute_offsets
; jm02feb7uofa

; mpage -l -1 -W145 offset_table_output.txt > offset_table_output.ps
    
    readcol, 'offset_table_input.txt', galaxy, dec, pa, gsize, offset, $
      scanlen20, scantime20, exp20, scanlen55, scantime20, exp55, $
      format='A,A,D,D,D,D,L,D,D,L,D', skip=1, /silent
    nobjects = n_elements(galaxy)
    
;   gsize = gsize*60.0          ; galaxy size [arcsec]
;   ssize = 2.0*130.0*15.0/18.0 ; approximate slit size [arcsec]

; 55" scans
    
; shift East then West

    dinfo55_E = driftscan(pa,dec,scanlen55,exp55,scantime=scantime20,slit_offset=offset)
    dinfo55_W = driftscan(pa,dec,scanlen55,exp55,scantime=scantime55,slit_offset=-offset)

; 20" scans of the center

    dinfo20 = driftscan(pa,dec,scanlen20,exp20,scantime=scantime20)
    
; write out 

    openw, 100, 'offset_table_output.txt'
    printf, 100, '# Note:  The PA increases from NORTH to EAST.  EAST is positive.'
    printf, 100, '# '+im_today()
    printf, 100, '# EAST scan'
    
; EAST
    
    outinfo_E = {galaxy: '', pa: 0.0D, size: 0.0D, offset: 0.0D, ra_rate_55: 0.0D, $
                 dec_rate_55: 0.0D, ra_offset_55: 0.0D, dec_offset_55: 0.0D}
    outinfo_E = replicate(outinfo_E,nobjects)

    outinfo_E.galaxy = galaxy
    outinfo_E.pa = pa
    outinfo_E.size = gsize
    outinfo_E.offset = offset

    outinfo_E.ra_rate_55 = dinfo55_E.drift_ra
    outinfo_E.dec_rate_55 = dinfo55_E.drift_dec
    outinfo_E.ra_offset_55 = dinfo55_E.offset_ra
    outinfo_E.dec_offset_55 = dinfo55_E.offset_dec

    print_struct, outinfo_E, lun_out=100

    printf, 100, ' '
    printf, 100, '# WEST scan'

; WEST

    outinfo_W = {galaxy: '', pa: 0.0D, size: 0.0D, offset: 0.0D, ra_rate_55: 0.0D, $
                 dec_rate_55: 0.0D, ra_offset_55: 0.0D, dec_offset_55: 0.0D}
    outinfo_W = replicate(outinfo_W,nobjects)

    outinfo_W.galaxy = galaxy
    outinfo_W.pa = pa
    outinfo_W.size = gsize
    outinfo_W.offset = -offset

    outinfo_W.ra_rate_55 = dinfo55_W.drift_ra
    outinfo_W.dec_rate_55 = dinfo55_W.drift_dec
    outinfo_W.ra_offset_55 = dinfo55_W.offset_ra
    outinfo_W.dec_offset_55 = dinfo55_W.offset_dec

    print_struct, outinfo_W, lun_out=100

; 20" scans

    printf, 100, ' '
    printf, 100, '# 20" Central Scans'

    outinfo = {galaxy: '', pa: 0.0D, ra_rate_20: 0.0D, $
               dec_rate_20: 0.0D, ra_offset_20: 0.0D, dec_offset_20: 0.0D}
    outinfo = replicate(outinfo,nobjects)
    outinfo.galaxy = galaxy
    outinfo.pa = pa

    outinfo.ra_rate_20 = dinfo20.drift_ra
    outinfo.dec_rate_20 = dinfo20.drift_dec
    outinfo.ra_offset_20 = dinfo20.offset_ra
    outinfo.dec_offset_20 = dinfo20.offset_dec

    print_struct, outinfo, lun_out=100
    print_struct, outinfo
    
    close, 100
    
    spawn, ['mpage -l -1 -W145 offset_table_output.txt > offset_table_output.ps']

return
end
