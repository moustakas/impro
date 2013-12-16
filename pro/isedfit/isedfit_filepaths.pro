;+
; NAME:
;   ISEDFIT_FILEPATHS()
; PURPOSE:
;   Given a parameters structure (see READ_ISEDFIT_PARAMFILE), pack
;   all the requisite paths and filenames into a structure.  This
;   routine is a support routine and in general should not be called
;   on its own.
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 11, NYU
;   jm13aug09siena - updated to the latest data model 
;-

function isedfit_filepaths, params, isedfit_dir=isedfit_dir, $
  montegrids_dir=montegrids_dir, outprefix=outprefix1, $
  sed_pdffile=sed_pdffile1, zcolor_pdffile=zcolor_pdffile1, $
  colorcolor_pdffile=colorcolor_pdffile1, priors_pdffile=priors_pdffile, $
  photoz_pdffile=photoz_pdffile1, band_shift=band_shift

    if (n_elements(params) eq 0) then $
      message, 'PARAMS input required'
    if (n_elements(params) ne 1) then $
      message, 'PARAMS must be a scalar'

    if n_elements(isedfit_dir) eq 0 then isedfit_dir = get_pwd()
    if n_elements(montegrids_dir) eq 0 then montegrids_dir = get_pwd()+'montegrids/'

; some handy variables    
    prefix = strtrim(params.prefix,2)
    if tag_exist(params,'MAXOLD') then if params.maxold then $
      suffix = '_maxold' else suffix = '' else suffix = ''

    if (n_elements(outprefix1) eq 0) then thisprefix = prefix else $
      thisprefix = outprefix1

    spsmodels = strtrim(params.spsmodels,2)
    imf = strtrim(params.imf,2)
    redcurve = strtrim(params.redcurve,2)
    sfhgridstring = 'sfhgrid'+string(params.sfhgrid,format='(I2.2)')

; montegrids output directory and output files
    montegrids_fullpath = montegrids_dir+sfhgridstring+'/'+spsmodels+'/'+redcurve+'/'
    montegrids_montefile = montegrids_fullpath+prefix+'_'+$
      spsmodels+'_'+imf+'_montegrid.fits'
    montegrids_chunkfiles = montegrids_fullpath+prefix+'_'+spsmodels+'_'+imf+$
      '_chunk_'+string(lindgen(params.nmodelchunk)+1,format='(I4.4)')+'.fits'

; convolved photometry (ISEDFIT_MODELS) output directory and files 
    models_fullpath = isedfit_dir+sfhgridstring+'/'+redcurve+'/'
    models_chunkfiles = models_fullpath+prefix+'_'+spsmodels+'_'+imf+$
      '_chunk_'+string(lindgen(params.nmodelchunk)+1,format='(I4.4)')+'.fits'
    
; iSEDfit output file names    
    outfile = thisprefix+'_'+spsmodels+'_'+imf+'_'+redcurve+'_'+sfhgridstring+suffix
    postfile = outfile+'_post'

    outfile_photoz = thisprefix+'_photoz_'+spsmodels+'_'+imf+'_'+redcurve+'_'+sfhgridstring+suffix
    postfile_photoz = outfile_photoz+'_post'

    if n_elements(band_shift) eq 0 then band_shift = 0.0
    zbandshift = 'z'+string(band_shift,format='(F3.1)')
    kcorrfile = outfile+'_kcorr.'+zbandshift
    
; QAplot file names    
    if n_elements(sed_pdffile1) eq 0 then begin
       sed_psfile = 'qaplot_sed_'+outfile+'.ps'
       sed_pdffile = repstr(sed_psfile,'.ps','.pdf')
    endif else begin
       sed_pdffile = file_basename(sed_pdffile1)
       sed_psfile = repstr(sed_pdffile,'.pdf','.ps')
    endelse

    if n_elements(zcolor_pdffile1) eq 0 then begin
       zcolor_psfile = 'qaplot_zcolor_'+outfile+'.ps'
       zcolor_pdffile = repstr(zcolor_psfile,'.ps','.pdf')
    endif else begin
       zcolor_pdffile = file_basename(zcolor_pdffile1)
       zcolor_psfile = repstr(zcolor_pdffile,'.pdf','.ps')
    endelse

    if n_elements(colorcolor_pdffile1) eq 0 then begin
       colorcolor_psfile = 'qaplot_colorcolor_'+outfile+'.ps'
       colorcolor_pdffile = repstr(colorcolor_psfile,'.ps','.pdf')
    endif else begin
       colorcolor_pdffile = file_basename(colorcolor_pdffile1)
       colorcolor_psfile = repstr(colorcolor_pdffile,'.pdf','.ps')
    endelse

    if n_elements(priors_pdffile1) eq 0 then begin
       priors_psfile = 'qaplot_priors_'+outfile+'.ps'
       priors_pdffile = repstr(priors_psfile,'.ps','.pdf')
    endif else begin
       priors_pdffile = file_basename(priors_pdffile1)
       priors_psfile = repstr(priors_pdffile,'.pdf','.ps')
    endelse

    if n_elements(photoz_pdffile1) eq 0 then begin
       photoz_psfile = 'qaplot_photoz_'+outfile+'.ps'
       photoz_pdffile = repstr(photoz_psfile,'.ps','.pdf')
    endif else begin
       photoz_pdffile = file_basename(photoz_pdffile1)
       photoz_psfile = repstr(photoz_pdffile,'.pdf','.ps')
    endelse

; file paths and filenames
    filepaths = {$
      isedfit_dir:               isedfit_dir,        $
      models_fullpath:           models_fullpath,         $
      models_chunkfiles:         models_chunkfiles,  $
      isedfit_outfile:           outfile+'.fits',  $
      isedfit_photoz_outfile:    outfile_photoz+'.fits',  $
      post_outfile:              postfile+'.fits',  $
      post_photoz_outfile:       postfile_photoz+'.fits',  $
      kcorr_outfile:             kcorrfile+'.fits',  $
      qaplot_sed_psfile:         sed_psfile,$
      qaplot_sed_pdffile:        sed_pdffile,$
      qaplot_zcolor_psfile:      zcolor_psfile,$
      qaplot_zcolor_pdffile:     zcolor_pdffile,$
      qaplot_colorcolor_psfile:  colorcolor_psfile,$
      qaplot_colorcolor_pdffile: colorcolor_pdffile,$
      qaplot_priors_psfile:      priors_psfile,$
      qaplot_priors_pdffile:     priors_pdffile,$
      qaplot_photoz_psfile:      photoz_psfile,$
      qaplot_photoz_pdffile:     photoz_pdffile,$
      montegrids_fullpath:       montegrids_fullpath,$
      montegrids_montefile:      montegrids_montefile,$
      montegrids_chunkfiles:     montegrids_chunkfiles}

return, filepaths
end
