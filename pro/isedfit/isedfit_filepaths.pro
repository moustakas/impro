;+
; NAME:
;   ISEDFIT_FILEPATHS()
;
; PURPOSE:
;   Given a parameters structure (see READ_ISEDFIT_PARAMFILE), pack
;   all the requisite paths and filenames into a structure.
;
; INPUTS: 
;   params - 
;
; OPTIONAL INPUTS: 
;   outprefix - 
;   isedpath - 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 11, NYU 
;
; Copyright (C) 2009, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function isedfit_filepaths, params, super=super, outprefix=outprefix1, $
  isedpath=isedpath, ngalaxy=ngalaxy, ngalchunk=ngalchunk, $
  galchunksize=galchunksize, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir

    if (n_elements(isedpath) eq 0L) then isedpath = './'
    if (n_elements(params) eq 0L) then $
      message, 'PARAMS input required'
    if (n_elements(galchunksize) eq 0L) then $
      galchunksize = 500L ; can play with this number

    if (n_elements(isedfit_sfhgrid_dir) eq 0) then isedfit_sfhgrid_dir = $
      '${ISEDFIT_SFHGRID_DIR}/'
    
; read the output from build_isedfit_sfhgrid based on the specified
; spectral synthesis models and SFH grid
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    redcurvestring = strtrim(redcurve2string(super.redcurve),2)
    sfhgrid = super.sfhgrid

;   synthmodels = strtrim(params.synthmodels,2)
;   imf = strtrim(params.imf,2)
;   redcurvestring = strtrim(params.redcurve,2)
;   sfhgrid = params.sfhgrid
    if (n_elements(sfhgrid) ne 1) then $
      message, 'SFHGRID must be a scalar'

; if CHUNKINFO file does not exist, try dropping the reddening curve 
    sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')
    sfhgridpath = isedfit_sfhgrid_dir+sfhgridstring+$
      '/'+synthmodels+'/'+redcurvestring+'/'
    chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits.gz'
    if (file_test(chunkinfofile) eq 0) then begin
       redcurvestring = ''
       sfhgridpath = isedfit_sfhgrid_dir+sfhgridstring+$
         '/'+synthmodels+'/'
       chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits.gz'
    endif
    if (file_test(chunkinfofile) eq 0) then $
      message, 'CHUNKINFOFILE '+chunkinfofile+' not found!'

    modelspath = isedpath+sfhgridstring+'/'+redcurvestring+'/'

; needs error checking here to make sure the directories exist    

; output file names prefix; PREFIX overrides everything; include an
; additional suffix if MAXOLD=1
    if (tag_exist(params,'PREFIX') eq 0) then prefix = 'isedfit' else $
      prefix = strtrim(params.prefix,2)
    if tag_exist(params,'MAXOLD') then if params.maxold then $
      suffix = '_maxold' else suffix = '' else suffix = ''
    
    if (n_elements(outprefix1) eq 0) then thisprefix = prefix else $
      thisprefix = outprefix1
    if (redcurvestring eq '') then $
      outfile = thisprefix+'_'+synthmodels+'_'+imf+'_'+sfhgridstring+suffix else $
      outfile = thisprefix+'_'+synthmodels+'_'+imf+'_'+redcurvestring+'_'+sfhgridstring+suffix
    postfile = outfile+'_post'
    kcorrfile = outfile+'_kcorr'
    psfile = 'qaplot_'+outfile+'.ps'
    chi2grid_psfile = 'qaplot_chi2grid_'+outfile+'.ps'

    modelsprefix = prefix+'_'+synthmodels
;   modelsprefix = prefix+'_'+imf+'_'+sfhgridstring+'_models'

; the galaxy sample must be split into chunks because the arrays can
; be pretty memory intensive; unfortunately, the following code has to
; go here (including the dependence on NGAL), instead of in ISEDFIT,
; because of the dependence on the number of output chi2grid files;
; this output is also needed by ISEDFIT_CHI2GRID_QAPLOT
    if (n_elements(ngalaxy) ne 0) then begin
       ngalchunk = ceil(ngalaxy/float(galchunksize))
; full chi2 grid file names
       len = strtrim(strlen(strtrim(ngalchunk,2))>3,2)
       chi2grid_gchunkfiles = modelsprefix+imf+'_chi2grid_gchunk_'+$
         string(lindgen(ngalchunk)+1,format='(I'+len+'.'+len+')')+'.fits'
;      modelgrid_gchunkfiles = modelsprefix+imf+'_modelgrid_gchunk_'+$
;        string(lindgen(ngalchunk)+1,format='(I'+len+'.'+len+')')+'.fits'
    endif else begin
       chi2grid_gchunkfiles = ''
       modelgrid_gchunkfiles = ''
    endelse

; file paths and filenames
    if (file_test(chunkinfofile,/regular) eq 0L) then $
      message, 'No SFH grid information file found!'
    chunkinfo = mrdfits(chunkinfofile,1,/silent)
    sfhgrid_chunkfiles = sfhgridpath+strtrim(chunkinfo.chunkfiles,2)

    sfhgrid_chunksuffix = repstr(file_basename(sfhgrid_chunkfiles),'.fits.gz','')
    isedfit_models_chunkfiles = modelsprefix+'_'+sfhgrid_chunksuffix+'.fits'

    filepaths = {$
      synthmodels:               synthmodels,        $
      sfhgrid:                   sfhgrid,            $
      isedpath:                    isedpath,             $
      modelspath:                modelspath,         $
      sfhgrid_dir:               sfhgridpath,        $
      sfhgrid_chunkinfo:         chunkinfofile,      $ 
      sfhgrid_chunkfiles:        sfhgrid_chunkfiles, $
      chi2grid_gchunkfiles:      chi2grid_gchunkfiles,$
      isedfit_models_chunkfiles: isedfit_models_chunkfiles,  $
      isedfit_outfile:           outfile+'.fits',  $
      post_outfile:              postfile+'.fits',  $
      kcorr_outfile:             kcorrfile+'.fits',  $
      qaplot_psfile:             psfile,$
      qaplot_chi2grid_psfile:    chi2grid_psfile}

return, filepaths
end
