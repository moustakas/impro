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
  pdffile=pdffile1, band_shift=band_shift

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

    synthmodels = strtrim(params.synthmodels,2)
    imf = strtrim(params.imf,2)
    redcurve = strtrim(params.redcurve,2)
    sfhgridstring = 'sfhgrid'+string(params.sfhgrid,format='(I2.2)')

; montegrids output directory and output files; only
; ISEDFIT_MONTEGRIDS, ISEDFIT_MODELS, ISEDFIT_RESTORE, and possibly
; ISEDFIT_RECONSTRUCT_POSTERIOR should need the CHUNKINFO file,
; so don't crash if we don't find it (just move on)
    montegrids_fullpath = montegrids_dir+sfhgridstring+'/'+synthmodels+'/'+redcurve+'/'
    montegrids_chunkinfofile = montegrids_fullpath+imf+'_chunkinfo.fits.gz'
    if file_test(montegrids_chunkinfofile) eq 0 then begin
;      message, 'CHUNKINFOFILE '+montegrids_chunkinfofile+' not found!'
       montegrids_chunkfiles = ''
    endif else begin
       chunkinfo = mrdfits(montegrids_chunkinfofile,1,/silent)
       montegrids_chunkfiles = montegrids_fullpath+strtrim(chunkinfo.chunkfiles,2)
    endelse

; convolved photometry (ISEDFIT_MODELS) output directory; these files
; are missing when ISEDFIT_MODELS is first called, at which point
; MONTEGRID_DIR is set and we can construct the names from
; MONTEGRIDS_CHUNKFILES 
    models_fullpath = isedfit_dir+sfhgridstring+'/'+redcurve+'/'
    models_chunkfiles = file_search(models_fullpath+prefix+'_'+synthmodels+'_'+$
      imf+'_chunk_????.fits.gz',count=nfile)
    if (nfile eq 0) then models_chunkfiles = models_fullpath+thisprefix+'_'+$
      synthmodels+'_'+file_basename(montegrids_chunkfiles)
    models_chunkfiles = repstr(models_chunkfiles,'.gz','')
    
; iSEDfit output file names    
    outfile = thisprefix+'_'+synthmodels+'_'+imf+'_'+redcurve+'_'+sfhgridstring+suffix
    postfile = outfile+'_post'
    if n_elements(pdffile1) eq 0 then begin
       psfile = 'qaplot_'+outfile+'.ps'
       pdffile = repstr(psfile,'.ps','.pdf')
    endif else begin
       pdffile = file_basename(pdffile1)
       psfile = repstr(pdffile,'.pdf','.ps')
    endelse

    if n_elements(band_shift) eq 0 then band_shift = 0.0
    zbandshift = 'z'+string(band_shift,format='(F3.1)')
    kcorrfile = outfile+'_kcorr.'+zbandshift
    
; file paths and filenames
    filepaths = {$
      isedfit_dir:              isedfit_dir,        $
      models_fullpath:          models_fullpath,         $
      models_chunkfiles:        models_chunkfiles,  $
      isedfit_outfile:          outfile+'.fits',  $
      post_outfile:             postfile+'.fits',  $
      kcorr_outfile:            kcorrfile+'.fits',  $
      qaplot_psfile:            psfile,$
      qaplot_pdffile:           pdffile,$
      montegrids_fullpath:      montegrids_fullpath,$
      montegrids_chunkinfofile: montegrids_chunkinfofile,$ 
      montegrids_chunkfiles:    montegrids_chunkfiles}

return, filepaths
end
