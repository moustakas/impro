function im_read_ssp, pegfile, bc03=bc03, pegase=pegase, _extra=extra
; jm09feb07nyu - wrapper on reading in various SSPs in a standard
; format

    lsun = 3.826D33         ; [erg/s]
    dist = 10.0*3.085678D18 ; fiducial distance [10 pc in cm]

    ssp = -1
    
    if keyword_set(bc03) then begin
       ssp = im_read_bc03(_extra=extra)
       ssp.flux = lsun*ssp.flux/(4.0*!dpi*dist^2.0) ; [erg/s/cm2/A/M_sun]
    endif
    
    if keyword_set(pegase) then begin
       ssp = im_read_peg(pegfile,_extra=extra)
    endif

return, ssp
end
    
