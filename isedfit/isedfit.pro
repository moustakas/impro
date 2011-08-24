;+
; NAME:
;   ISEDFIT
;
; PURPOSE:
;   Given photometry in NFILT bands, and the output from
;   ISEDFIT_MODELS, compute the best-fitting model parameters
;   (mass, reddening, age, tau value, and stellar metallicity). 
;
; INPUTS:
;   maggies     - photometry [NFILT,NGALAXY]
;   ivarmaggies - inverse variance array for MAGGIES [NFILT,NGALAXY]  
;   zobj        - galaxy redshift [NGALAXY] 
;
; OPTIONAL INPUTS:
;   filterlist    - use photometry only in this list of filters;
;                   must be all or a subset of the filters used in 
;                   ISEDFIT_MODELS
;   datapath      - must match ISEDFIT_MODELS [default CWD()]
;   modelsprefix  - see ISEDFIT_MODELS (default 'isedfit_models')
;   isedfitprefix - output prefix; may be different than
;                   MODELSPREFIX if using only a subset of the 
;                   filters (see FILTERLIST) (default 'isedfit') 
;   galaxy        - optional galaxy name [NGALAXY] 
;   nmonte        - number of Monte Carlo realizations [NMONTE] 
;   nminphot      - require good photometry in at least NMINPHOT
;                   bandpasses to perform the fitting
;
; KEYWORD PARAMETERS:
;   milkyway   - use the O'Donnell (1994) Milky Way extinction
;                curve (default is to use Charlot & Fall 2000)
;   smc        - use the Gordon et al. 2003 SMC bar extinction
;                curve (default is to use Charlot & Fall 2000)
;   calzetti   - use the Calzetti 2000 continuum attenuation
;                curve for starburst galaxies (default is to use
;                Charlot & Fall 2000)
;   maxold     - include a second, maximally old component
;   noagelimit - allow solutions with ages that are older than the
;                age of the universe at the redshift of the object
;   write      - write out
;
; OUTPUTS:
;   result      - output data structure
;   result_info - corresponding information structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   No error checking is done to verify that the MAGGIES through
;   the appropriate filters were computed in ISEDFIT_MODELS,
;   relative to the MAGGIES array that is given to this routine.
;   Garbage in, garbage out.
;
; TODO:
;   If a subset of the available filters are used, then still
;   store the best-fitting magnitudes of the filters not used so
;   they can be plotted in ISEDFIT_QAPLOT.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Feb 12, U of A
;   jm05feb22uofa - better memory management
;   jm05mar22uofa - added optional input FILTERLIST, which can be
;                   a subset of the filters used in ISEDFIT_MODELS 
;   jm05may24uofa - bug fix in identifying multiple chi2 minima
;   jm05aug04uofa - write out RESULT and CHI2RESULT after every
;                   iteration in case the code crashes before
;                   finishing 
;   jm06mar02uofa - major re-write & upgrade
;   jm06mar15uofa - additional error checking for infinite values
;                   in the input data, and verifying that the
;                   filter dimensions of MAGGIES matches
;                   FILTERLIST
;   jm06mar21uofa - added MAXOLD and NOAGELIMIT keywords 
;   jm06jul13uofa - various bug fixes
;   jm07mar15nyu  - further developments
;
; Copyright (C) 2005-2007, John Moustakas
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
      ebv:            -1.0,$
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
      ebv_avg:      -1.0,$
      mu_avg:       -1.0,$

      mass_err:     -1.0,$
      age_err:      -1.0,$
      sfr_err:      -1.0,$
      sfr100_err:   -1.0,$
      b100_err:     -1.0,$
      tau_err:      -1.0,$
      Z_err:        -1.0,$
      ebv_err:      -1.0,$
      mu_err:       -1.0,$

      mass_50:     -1.0,$
      age_50:      -1.0,$
      sfr_50:      -1.0,$ ; instantaneous
      sfr100_50:   -1.0,$ ; 100 Myr
      b100_50:     -1.0,$
      tau_50:      -1.0,$
      Z_50:        -1.0,$
      ebv_50:      -1.0,$
      mu_50:       -1.0,$

      mass_eff_err:   -1.0,$
      age_eff_err:    -1.0,$
      sfr_eff_err:    -1.0,$
      sfr100_eff_err: -1.0,$
      b100_eff_err:   -1.0,$
      tau_eff_err:    -1.0,$
      Z_eff_err:      -1.0,$
      ebv_eff_err:    -1.0,$
      mu_eff_err:     -1.0}

;     mass_mode:     -1.0,$
;     age_mode:      -1.0,$
;     sfr_mode:      -1.0,$ ; instantaneous
;     sfr100_mode:   -1.0,$ ; 100 Myr
;     b100_mode:     -1.0,$
;     tau_mode:      -1.0,$
;     Z_mode:        -1.0,$
;     ebv_mode:      -1.0,$
;     mu_mode:       -1.0}

    isedfit = struct_addtags(struct_addtags(best,qmed),isedfit1)
    isedfit = replicate(temporary(isedfit),ngal)

; initialize the posterior distribution structure
    isedfit_post = {$
      draws:     lonarr(ndraw)-1,$
      scale:     fltarr(ndraw)-1,$
      scale_err: fltarr(ndraw)-1}
;     mass:   fltarr(ndraw)-1}
;     Z:      fltarr(ndraw),$
;     age:    fltarr(ndraw),$
;     tau:    fltarr(ndraw),$
;     ebv:    fltarr(ndraw),$
;     mu:     fltarr(ndraw),$
;     sfr:    fltarr(ndraw),$
;     sfr100: fltarr(ndraw),$
;     b100:   fltarr(ndraw)}
    isedfit_post = replicate(temporary(isedfit_post),ngal)
    
return, isedfit
end

pro isedfit, paramfile, maggies, ivarmaggies, zobj, isedfit, isedfit_post=isedfit_post, $
  params=params, iopath=iopath, nminphot=nminphot, galchunksize=galchunksize, $
  outprefix=outprefix, sfhgrid_paramfile=sfhgrid_paramfile, sfhgrid_basedir=sfhgrid_basedir, $
  index=index, allages=allages, clobber=clobber, debug=debug, write_chi2grid=write_chi2grid, $
  nowrite=nowrite, silent=silent

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

; SFHGRID and REDCURVE can be vectors; however, if SFHGRID=3 (no
; reddening) then ignore REDCURVE    
    nsfhgrid = n_elements(params.sfhgrid)
    nredcurve = n_elements(params.redcurve)
    if (nsfhgrid gt 1) or (nredcurve gt 1) then begin
       for ii = 0, nsfhgrid-1 do begin
          newparams1 = struct_trimtags(params,except='sfhgrid')
          newparams1 = struct_addtags(newparams1,{sfhgrid: params.sfhgrid[ii]})
          if (newparams1.sfhgrid eq 3) then nredcurve = 1
          for jj = 0, nredcurve-1 do begin
             newparams2 = struct_trimtags(newparams1,except='redcurve')
             newparams2 = struct_addtags(newparams2,{redcurve: params.redcurve[jj]})
             isedfit, paramfile, maggies, ivarmaggies, zobj, isedfit, isedfit_post=isedfit_post, $
               params=newparams2, iopath=iopath, nminphot=nminphot, $
               galchunksize=galchunksize, outprefix=outprefix, index=index, $
               sfhgrid_paramfile=sfhgrid_paramfile, sfhgrid_basedir=sfhgrid_basedir, $
               allages=allages, clobber=clobber, debug=debug, write_chi2grid=write_chi2grid, $
               nowrite=nowrite, silent=silent
          endfor
       endfor
       return
    endif

    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      ngalaxy=ngal,ngalchunk=ngalchunk,galchunksize=galchunksize,$
      sfhgrid_basedir=sfhgrid_basedir)

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
         debug=debug, write_chi2grid=write_chi2grid, /nowrite, silent=silent, $
         sfhgrid_paramfile=sfhgrid_paramfile, sfhgrid_basedir=sfhgrid_basedir
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
    if (min(zobj) lt min(redshift)) or $
      (max(zobj) gt max(redshift)) then begin
       splog, 'Need to rebuild model grids using a wider redshift grid!'
       return
    endif
    zindx = findex(redshift,zobj) ; used for interpolation

; do not allow the galaxy to be older than the age of the universe at
; ZOBJ starting from a maximum formation redshift z=10 [Gyr]
    maxage = getage(zobj,/gyr)-getage(10.0,/gyr) 

; initialize the output structure(s)
    isedfit = init_isedfit(ngal,nfilt,params.sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile,$
      isedfit_post=isedfit_post)
    isedfit.isedfit_id = lindgen(ngal)
    isedfit.maggies = maggies
    isedfit.ivarmaggies = ivarmaggies
    isedfit.zobj = zobj

; loop on each galaxy chunk    
    t1 = systime(1)
;   for gchunk = 2, ngalchunk-1L do begin
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
         modelgrid,fullgrid,isedfit_post=temp_isedfit_post,debug=debug)
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
