pro make_bok_catalog, ra, de, mag, object, epoch, datapath=datapath, catname=catname
; jm02jun24uofa
; generate a source catalog that is compatible with the 90"

    nobject = n_elements(ra)
    id = lindgen(nobject)
    
    if n_elements(datapath) eq 0L then begin
       spawn, ['pwd'], datapath
       datapath = datapath[0]+'/'
    endif

    if n_elements(catname) eq 0L then catname = 'catalog'
    
    dtor = !dpi/180.0D ; degrees --> radians
    
; convert RA and DEC to radians

    rra = 15D0*dtor*im_hms2dec(ra) ; [rad]
    rde = dtor*im_hms2dec(de)      ; [rad]

; define the output formats

    frmtpos = '(x,I3,x,I4,x,F12.10,x,A1,F12.10,6x,A1,6x,A1,x,A9,56x,F7.2)'
    frmtneg = '(x,I3,x,I4,x,F12.10,x,F13.10,6x,A1,6x,A1,x,A9,56x,F7.2)'

    rdepos = where(rde gt 0.0,nrdepos,comp=rdeneg,ncomp=nrdeneg)
    
    openw, lun, datapath+catname+'.cat', /get_lun
    for i = 0L, nobject-1L do begin
       if rde[i] gt 0.0 then printf, lun, i, mag[i], rra[i], '+', rde[i], '0', '0', object[i], epoch[i], format=frmtpos else $
         printf, lun, i, mag[i], rra[i], rde[i], '0', '0', object[i], epoch[i], format=frmtneg
    endfor
    free_lun, lun

return
end
