;+
; NAME:
;   ISEDFIT
;
; PURPOSE:
;   Model the spectral energy distributions of galaxies using
;   population synthesis models to infer their physical properties. 
;
; INPUTS:
;   paramfile - iSEDfit parameter file see README_ISEDFIT for details) 
;   maggies - photometry [NFILT,NGAL]
;   ivarmaggies - inverse variance array for MAGGIES [NFILT,NGAL]  
;   zobj - galaxy redshift [NGAL] 
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     PARAMFILE (over-rides PARAMFILE)
;   iopath - full path name to the input and output files (default ./) 
;   nminphot - require at least NMINPHOT bandpasses of well-measured
;     photometry (i.e., excluding upper limits) before fitting
;     (default 3)  
;   galchunksize - split the sample into GALCHUNKSIZE sized chunks,
;     which is necessary if the sample is very large (default 5000)
;   outprefix - optionally write out files with a different prefix
;     from that specified in PARAMFILE (or PARAMS)
;   sfhgrid_paramfile - parameter file name giving details on all the
;     available SFH grids (needed to initialize the output data
;     structure) 
;   isedfit_sfhgrid_dir - full pathname to the precomputed SFH grids 
;   index - use this optional input to fit a zero-indexed subset of
;     the full sample (default is to fit everything)
;
; KEYWORD PARAMETERS:
;   allages - allow solutions with ages that are older than the age of
;     the universe at the redshift of the object 
;   write_chi2grid - write out the full chi^2 grid across all the
;     models (these can be big files!)
;   silent - suppress messages to STDOUT
;   nowrite - do not write out any of the output files (generally not
;     recommended!) 
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;
; OUTPUTS:
;   
;
; OPTIONAL OUTPUTS:
;   isedfit - output data structure containing all the results; see
;     README_ISEDFIT for a detailed breakdown and explanation of all
;     the outputs
;   isedfit_post - output data structure containing the random draws
;     from the posterior distribution function, which can be used to
;     rebuild the posterior distributions of any of the output
;     parameters 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Sep 01, UCSD - I began writing iSEDfit in 2005
;     while at the U of A, adding updates off-and-on through 2007;
;     however, the code has evolved so much that the old modification
;     history became obsolete!  Future changes to the officially
;     released code will be documented here.
;
; Copyright (C) 2011, John Moustakas
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

function init_isedfit, ngal, nfilt, sfhgrid, sfhgrid_paramfile=sfhgrid_paramfile, $
  isedfit_post=isedfit_post
; ISEDFIT support routine - initialize the output structure 

    params = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ndraw = isedfit_ndraw() ; number of random draws

; compute the maximum number of bursts, if any
    if (params.pburst le 0.0) then nmaxburst = 0 else $
      nmaxburst = ceil((params.maxage-params.minage)/params.pburstinterval)

    burstarray1 = -1.0
    burstarray2 = -1.0
    if (nmaxburst gt 1) then begin
       burstarray1 = fltarr(nmaxburst)-1.0
       burstarray2 = dblarr(nmaxburst)-1D
    endif

    isedfit1 = {$
      isedfit_id:      -1L,$    ; unique ID number
      zobj:           -1.0,$    ; redshift
      maggies:     fltarr(nfilt),$ ; observed maggies
      ivarmaggies: fltarr(nfilt),$ ; corresponding inverse variance
      bestmaggies: fltarr(nfilt)}  ; best-fitting model photometry

; best-fit values (at the chi2 minimum); see BUILD_ISEDFIT_SFHGRID
    best = {$
      imf:              '',$
      chunkindx:        -1,$
      modelindx:       -1L,$
      ageindx:          -1,$

      tau:            -1.0,$
      Z:              -1.0,$
      av:             -1.0,$
      mu:             -1.0,$
      nburst:            0,$
      tauburst:        -1D,$ ; burst truncation time scale
      tburst:    burstarray2,$
      dtburst:   burstarray2,$
      fburst:    burstarray1,$

      scale:          -1.0,$ 
      scale_err:      -1.0,$ 
      mass:           -1.0,$ 
      age:            -1.0,$
      sfr:            -1.0,$ ; instantaneous
      sfr100:         -1.0,$ ; 100 Myr
      b100:           -1.0,$ ; averaged over the previous 100 Myr
      chi2:            1E6}  ; chi2 minimum

; median quantities and PDF quantiles
    qmed = {$
      mass_avg:     -1.0,$
      age_avg:      -1.0,$
      sfr_avg:      -1.0,$ ; instantaneous
      sfr100_avg:   -1.0,$ ; 100 Myr
      b100_avg:     -1.0,$
      tau_avg:      -1.0,$
      Z_avg:        -1.0,$
      av_avg:       -1.0,$
      mu_avg:       -1.0,$

      mass_50:     -1.0,$
      age_50:      -1.0,$
      sfr_50:      -1.0,$ ; instantaneous
      sfr100_50:   -1.0,$ ; 100 Myr
      b100_50:     -1.0,$
      tau_50:      -1.0,$
      Z_50:        -1.0,$
      av_50:       -1.0,$
      mu_50:       -1.0,$

      mass_err:     -1.0,$
      age_err:      -1.0,$
      sfr_err:      -1.0,$
      sfr100_err:   -1.0,$
      b100_err:     -1.0,$
      tau_err:      -1.0,$
      Z_err:        -1.0,$
      av_err:       -1.0,$
      mu_err:       -1.0}

;     mass_eff_err:   -1.0,$
;     age_eff_err:    -1.0,$
;     sfr_eff_err:    -1.0,$
;     sfr100_eff_err: -1.0,$
;     b100_eff_err:   -1.0,$
;     tau_eff_err:    -1.0,$
;     Z_eff_err:      -1.0,$
;     av_eff_err:     -1.0,$
;     mu_eff_err:     -1.0}

    isedfit = struct_addtags(struct_addtags(best,qmed),isedfit1)
    isedfit = replicate(temporary(isedfit),ngal)

; initialize the posterior distribution structure
    isedfit_post = {$
      draws:     lonarr(ndraw)-1,$
      scale:     fltarr(ndraw)-1,$
      scale_err: fltarr(ndraw)-1}
    isedfit_post = replicate(temporary(isedfit_post),ngal)
    
return, isedfit
end

pro isedfit, paramfile, maggies, ivarmaggies, zobj, isedfit, isedfit_post=isedfit_post, $
  params=params, iopath=iopath, nminphot=nminphot, galchunksize=galchunksize, $
  outprefix=outprefix, sfhgrid_paramfile=sfhgrid_paramfile, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
  index=index, allages=allages, write_chi2grid=write_chi2grid, silent=silent, $
  nowrite=nowrite, clobber=clobber

    if (n_elements(paramfile) eq 0) and $
      (n_elements(params) eq 0) then begin
       doc_library, 'isedfit'
       return
    endif

    ndim = size(maggies,/n_dim)
    dims = size(maggies,/dim)
    if (ndim eq 1) then ngal = 1 else ngal = dims[1]  ; number of galaxies
    nfilt = dims[0] ; number of filters

; error checking on the input photometry    
    nmaggies = n_elements(maggies)
    nivarmaggies = n_elements(ivarmaggies)
    nzobj = n_elements(zobj)
    
    if (nmaggies eq 0L) or (nivarmaggies eq 0L) or $
      (nzobj eq 0L) then begin
       doc_library, 'isedfit'
       return
    endif

    ndim = size(maggies,/n_dimension)
    if (ndim ne 2L) then begin ; reform into a 2D array
       maggies = reform(maggies,n_elements(maggies),1)
       ivarmaggies = reform(ivarmaggies,n_elements(maggies),1)
    endif

    if (n_elements(maggies) ne n_elements(ivarmaggies)) then begin
       splog, 'MAGGIES and IVARMAGGIES do not match'
       return
    endif
    if (nzobj ne ngal) then begin
       splog, 'MAGGIES and ZOBJ do not match'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(iopath) eq 0) then iopath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)

; SFHGRID can be a vector
    nsfhgrid = n_elements(params.sfhgrid)
    if (nsfhgrid gt 1) then begin
       for ii = 0, nsfhgrid-1 do begin
          newparams1 = struct_trimtags(params,except='sfhgrid')
          newparams1 = struct_addtags(newparams1,{sfhgrid: params.sfhgrid[ii]})
          isedfit, paramfile, maggies, ivarmaggies, zobj, isedfit, isedfit_post=isedfit_post, $
            params=newparams1, iopath=iopath, nminphot=nminphot, $
            galchunksize=galchunksize, outprefix=outprefix, index=index, $
            sfhgrid_paramfile=sfhgrid_paramfile, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            allages=allages, clobber=clobber, write_chi2grid=write_chi2grid, $
            nowrite=nowrite, silent=silent
       endfor
       return
    endif

    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      ngalaxy=ngal,ngalchunk=ngalchunk,galchunksize=galchunksize,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    outfile = fp.iopath+fp.isedfit_outfile
    if file_test(outfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) and $
      (keyword_set(nowrite) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif

; fit the requested subset of objects and return
    if (n_elements(index) ne 0L) then begin
       isedfit, paramfile, maggies[*,index], ivarmaggies[*,index], zobj[index], $
         isedfit1, isedfit_post=isedfit_post1, params=params, nminphot=nminphot, $
         outprefix=outprefix, allages=allages, iopath=iopath, clobber=clobber, $
         write_chi2grid=write_chi2grid, /nowrite, silent=silent, $
         sfhgrid_paramfile=sfhgrid_paramfile, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       isedfit = init_isedfit(ngal,nfilt,params.sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile,$
         isedfit_post=isedfit_post)
       isedfit[index] = isedfit1
       isedfit_post[index] = isedfit_post1
       if (keyword_set(nowrite) eq 0) then begin
          im_mwrfits, isedfit, outfile, /clobber
          im_mwrfits, isedfit_post, fp.iopath+fp.post_outfile, /clobber
       endif
       return
    endif

; additional error checking    
    if (total(finite(maggies) eq 0B) ne 0.0) or $
      (total(finite(ivarmaggies) eq 0B) ne 0.0) then begin
       splog, 'MAGGIES and IVARMAGGIES cannot have infinite values'
       return
    endif
    if (total(zobj le 0.0) ne 0.0) then begin
       splog, 'ZOBJ should all be positive'
       return
    endif
    if (n_elements(nminphot) eq 0L) then nminphot = 3

    if (keyword_set(silent) eq 0) then begin
       splog, 'SYNTHMODELS='+params.synthmodels+', '+$
         'IMF='+params.imf+', '+'SFHGRID='+$
         string(params.sfhgrid,format='(I2.2)')+', '+$
         'REDCURVE='+strtrim(params.redcurve,2)
    endif

; filters and redshift grid
    filterlist = strtrim(params.filterlist,2)
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filterlist)

    redshift = params.redshift
    if (min(zobj)-min(redshift) lt -1E-3) or $
      (max(zobj)-max(redshift) gt 1E-3) then begin
       splog, 'Need to rebuild model grids using a wider redshift grid!'
       return
    endif
    zindx = findex(redshift,zobj) ; used for interpolation

; do not allow the galaxy to be older than the age of the universe at
; ZOBJ starting from a maximum formation redshift z=10 [Gyr]
    maxage = getage(zobj,/gyr)
;   maxage = getage(zobj,/gyr)-getage(10.0,/gyr) 

; initialize the output structure(s)
    isedfit = init_isedfit(ngal,nfilt,params.sfhgrid,isedfit_post=isedfit_post,$
      sfhgrid_paramfile=sfhgrid_paramfile)
    isedfit.isedfit_id = lindgen(ngal)
    isedfit.maggies = maggies
    isedfit.ivarmaggies = ivarmaggies
    isedfit.zobj = zobj

; loop on each galaxy chunk    
    t1 = systime(1)
    for gchunk = 0L, ngalchunk-1L do begin
       g1 = gchunk*galchunksize
       g2 = ((gchunk*galchunksize+galchunksize)<ngal)-1L
       gnthese = g2-g1+1L
       gthese = lindgen(gnthese)+g1
; loop on each "chunk" of output from ISEDFIT_MODELS
       nchunk = n_elements(fp.isedfit_models_chunkfiles)
       t0 = systime(1)
       for ichunk = 0, nchunk-1 do begin
          chunkfile = fp.modelspath+fp.isedfit_models_chunkfiles[ichunk]+'.gz'
          if (keyword_set(silent) eq 0) then splog, 'Reading '+chunkfile
          chunkmodels = mrdfits(chunkfile,1,/silent)
          nmodel = n_elements(chunkmodels)
; compute chi2
          gridchunk = isedfit_compute_chi2(maggies[*,gthese],ivarmaggies[*,gthese],$
            chunkmodels,maxage[gthese],zindx[gthese],gchunk=gchunk,ngalchunk=ngalchunk,$
            ichunk=ichunk,nchunk=nchunk,nminphot=nminphot,allages=allages,silent=silent,$
            maxold=params.maxold)
          modelgrid1 = struct_trimtags(temporary(chunkmodels),except=['MODELMAGGIES'])
          if (ichunk eq 0) then begin
             fullgrid = temporary(gridchunk)
             modelgrid = temporary(modelgrid1)
;            nallmodel = temporary(nmodel)
          endif else begin
;            if ichunk+1 eq 6 then stop
             fullgrid = [[temporary(fullgrid)],[temporary(gridchunk)]]
             modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
;            nallmodel = nmodel+nallmodel
          endelse
       endfor ; model chunk
       if (keyword_set(silent) eq 0) then splog, format='("All ModelChunks = '+$
         '",G0," minutes          ")', (systime(1)-t0)/60.0
; optionally write out the full chi2 grid 
       if keyword_set(write_chi2grid) then begin
          im_mwrfits, fullgrid, clobber=clobber, $
            fp.modelspath+fp.chi2grid_gchunkfiles[gchunk]
;         im_mwrfits, modelgrid, clobber=clobber, $
;           fp.modelspath+fp.modelgrid_gchunkfiles[gchunk]
       endif
; minimize chi2
       if (keyword_set(silent) eq 0) then splog, 'Minimizing chi2...'
       t0 = systime(1)
       temp_isedfit_post = isedfit_post[gthese]
       isedfit[gthese] = isedfit_compute_posterior(isedfit[gthese],$
         modelgrid,fullgrid,isedfit_post=temp_isedfit_post)
       isedfit_post[gthese] = temporary(temp_isedfit_post) ; pass-by-value
       if (keyword_set(silent) eq 0) then splog, format='("Time = '+$
         '",G0," minutes")', (systime(1)-t0)/60.0
    endfor ; close GalaxyChunk
    if (keyword_set(silent) eq 0) then splog, format='("All GalaxyChunks = '+$
      '",G0," minutes        ")', (systime(1)-t1)/60.0

; write out the final structure and the full posterior distributions
    if (keyword_set(nowrite) eq 0) then begin
       im_mwrfits, isedfit, outfile, silent=silent, /clobber
       im_mwrfits, isedfit_post, fp.iopath+fp.post_outfile, silent=silent, /clobber
    endif

return
end
