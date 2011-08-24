function get_sfhgrid_colors, sfhgrid, synthmodels=synthmodels, imf=imf, $
  redcurve=redcurve, sfhgrid_paramfile=sfhgrid_paramfile, $
  isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, filterlist=filterlist, $
  band_shift=band_shift
; jm10jan28ucsd - take the output from BUILD_ISEDFIT_SFHGRID and
; synthesize photometry in an arbitrary set of bandpasses

    if (n_elements(sfhgrid) eq 0) then begin
       doc_library, 'get_sfhgrid_colors'
       return, -1
    endif

    nfilt = n_elements(filterlist)
    if (nfilt eq 0) then begin
       splog, 'FILTERLIST input required'
       return, -1
    endif

    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(imf) eq 0) then imf = 'chab' ; 'salp'
    if (n_elements(redcurve) eq 0) then redcurve = 0

; read the parameter file describing each of the grids and get the
; reddening curve
;   params = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    redcurvestring = redcurve2string(redcurve)

    if (n_elements(isedfit_sfhgrid_dir) eq 0) then isedfit_sfhgrid_dir = $
      '${ISEDFIT_SFHGRID_DIR}/'
    sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')
    sfhgridpath = isedfit_sfhgrid_dir+sfhgridstring+$
      '/'+synthmodels+'/'+redcurvestring+'/'

    chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits.gz'
    if (file_test(chunkinfofile) eq 0) then begin
       splog, 'No SFH grid information file found!'
       return, -1
    endif

; read the chunks and synthesize the photometry    
    chunkinfo = mrdfits(chunkinfofile,1,/silent)
    chunkfiles = sfhgridpath+strtrim(chunkinfo.chunkfiles,2)
    nchunk = n_elements(chunkfiles)

    for ichunk = 0L, nchunk-1L do begin
       splog, 'Reading '+chunkfiles[ichunk]
       info1 = mrdfits(chunkfiles[ichunk],1)
       nmodel = n_elements(info1)
       colors1 = struct_trimtags(info1,except=['wave','flux'])
       colors1 = struct_addtags(temporary(colors1),$
         replicate({abmag: fltarr(nfilt,n_elements(info1[0].age))},nmodel))
       
       kwave = k_lambda_to_edges(info1[0].wave)
       for ii = 0L, nmodel-1 do begin
          k_projection_table, rmatrix, info1[ii].flux, kwave, $
            /silent, zvals, filterlist, zmin=0.0, zmax=0.0, nz=1, $
            band_shift=band_shift
          colors1[ii].abmag = -2.5*alog10(transpose(reform(rmatrix)))
       endfor
       if (ichunk eq 0L) then colors = colors1 else $
         colors = [temporary(colors),colors1]
    endfor 
    
return, colors
end
