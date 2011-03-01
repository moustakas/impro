function read_sfhgrid_paramfile, sfhgrid, sfhgrid_paramfile=sfhgrid_paramfile
; jm10nov12ucsd - read the SFHGRID definitions parameter file

    if (n_elements(sfhgrid_paramfile) eq 0) then sfhgrid_paramfile = $
      getenv('IMPRO_DIR')+'/isedfit/isedfit_sfhgrid.par'
    if (file_test(sfhgrid_paramfile,/regular) eq 0L) then $
      message, 'PARAMFILE '+sfhgrid_paramfile+' not found'

    allparams = yanny_readone(sfhgrid_paramfile)

; build the parameter arrays for the specified SFHGRID
    if (n_elements(sfhgrid) ne 0) then begin
       match = where(allparams.sfhgrid eq sfhgrid)
       if (match[0] eq -1) then message, 'Please update '+sfhgrid_paramfile
       params = allparams[match]
       return, params
    endif

return, allparams
end
    
