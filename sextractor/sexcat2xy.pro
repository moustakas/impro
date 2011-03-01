pro sexcat2xy, root, datapath, catlist=catlist, outpath=outpath
; jm07jun26nyu - read in a SE catalog and output an xy region file for
;                DS9; the SE catalog must be in FITS_LDAC format

    if (n_elements(root) eq 0L) then root = ''
    if (n_elements(datapath) eq 0L) then datapath = './'
    if (n_elements(outpath) eq 0L) then outpath = datapath

    if (n_elements(catlist) eq 0L) then begin
       catlist = file_search(datapath+root+'*.cat',count=ncat)
       if (ncat eq 0L) then begin
          splog, 'No files found!'
          return
       endif
    endif else ncat = n_elements(catlist)

    for icat = 0L, ncat-1L do begin
       if (file_test(catlist[icat]) eq 0L) then begin
          splog, 'Catalog '+catlist[icat]+' not found'
          continue
       endif
       fits_info, catlist[icat], n_ext=next, /silent
       for iext = 0L, next/2L-1L do begin
          cat = mrdfits(catlist[icat],2*(iext+1L),/silent)
;         stars = lindgen(n_elements(cat))
;         stars = where(cat.class_star gt 0.8)
          stars = where(cat.flux_auto/cat.fluxerr_auto gt 10.0,nstars)
          if (nstars ne 0L) then begin
             reglist = outpath+file_basename(repstr(catlist[icat],'.cat','')+'.ext'+string(iext+1L,format='(I0)')+'.reg')
             splog, 'Writing '+reglist
             struct_print, struct_trimtags(cat[stars],select=['XWIN_IMAGE','YWIN_IMAGE']), file=reglist, /no_head
          endif else splog, 'No stars identified in catalog '+catlist[icat]+', extension '+string(iext+1L,format='(I0)')
       endfor 
    endfor 
    
return
end
    
