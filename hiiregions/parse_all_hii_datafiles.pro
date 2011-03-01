pro parse_all_hii_datafiles
; jm05may13uofa
; jm06jan18uofa - updated
; jm07dec17nyu - rewritten (generalized)

    path = hiiregions_path()

    spawn, 'find . -name "parse_*.pro" -print', all

    for iall = 0L, n_elements(all)-1L do begin
       splog, '##################################################'
       splog, 'Executing '+file_basename(all[iall])
       call_procedure, repstr(file_basename(all[iall]),'.pro','')
       splog, '##################################################'
    endfor

return
end
    
