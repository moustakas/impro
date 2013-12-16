;+
; NAME:
;   ISEDFIT_RECONSTRUCT_POSTERIOR()
;
; PURPOSE:
;   Reconstruct the posterior distribution(s) of the output
;   parameters.  This routine is used by READ_ISEDFIT() among other
;   routines but it can also be called independently. 
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the iSEDfit models and
;     output files should be written (default PWD=present working
;     directory) 
;   montegrids_dir - full directory path where the Monte Carlo grids
;     written by ISEDFIT_MONTEGRIDS can be found (default 'montegrids'
;     subdirectory of the PWD=present working directory)
;   index - use this optional input to restore a zero-indexed subset
;     of the full sample (default is to restore everything)
;   outprefix - optional output prefix string (see ISEDFIT) 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   result - data structure array containing lots of goodies 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Memory problems may be encountered if this routine is used to
;   reconstruct too many models/results.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 27, NYU - largely excised from
;     ISEDFIT_SED_QAPLOT and ISEDFIT_MEASURE
;   jm13aug09siena - updated to conform to the latest version of
;     iSEDfit 
;
; Copyright (C) 2007, 2013, John Moustakas
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

function isedfit_reconstruct_posterior, isedfit_paramfile, params=params, $
  thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
  index=index, outprefix=outprefix

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_reconstruct_posterior'
       return, -1
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if n_elements(params) eq 0 then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if n_elements(isedfit_dir) eq 0 then isedfit_dir = get_pwd()
    if n_elements(montegrids_dir) eq 0 then montegrids_dir = get_pwd()+'montegrids/'

; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          result1 = isedfit_reconstruct_posterior(params=params[ii],$
            isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
            index=index,outprefix=outprefix)
          if ii eq 0 then result = temporary(result1) else $
            result = [[temporary(result)],[temporary(result1)]]
       endfor
       return, result
    endif

    fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,outprefix=outprefix)

; restore the best-fitting results
    isedfile = fp.isedfit_dir+fp.isedfit_outfile+'.gz'
    if file_test(isedfile) eq 0 then begin
       splog, 'iSEDfit output file '+isedfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+isedfile
    isedfit1 = mrdfits(isedfile,1,/silent)

; only restore a subset of the objects    
    if (n_elements(index) eq 0L) then isedfit = temporary(isedfit1) else $
      isedfit = temporary(isedfit1[index])
    ngal = n_elements(isedfit)
    
; restore the posterior distribution file
    if file_test(fp.isedfit_dir+fp.post_outfile+'.gz') eq 0 then begin
       splog, 'POSTERIOR output file '+fp.isedfit_dir+fp.post_outfile+'.gz not found!'
       return, -1
    endif
    splog, 'Reading '+fp.isedfit_dir+fp.post_outfile+'.gz'
    post = mrdfits(fp.isedfit_dir+fp.post_outfile+'.gz',1,/silent,rows=index)
    ngal = n_elements(post)

; need all the chunks in memory!
    tags = ['chunkindx','modelindx','age','tau','zmetal','av','mu',$
      'mstar','sfr','sfr100','b100','b1000','sfrage','ewoii','ewoiiihb','ewniiha']
    nchunk = n_elements(fp.montegrids_chunkfiles)
    for jj = 0, nchunk-1 do begin
       chunkfile = fp.models_chunkfiles[jj]+'.gz'
       if file_test(chunkfile) eq 0 then message, 'Chunk file '+chunkfile+' not found!'
       modelgrid1 = mrdfits(chunkfile,1,/silent);,columns=tags)
       if (jj eq 0) then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor
    nmodel = n_elements(modelgrid)
    ndraw = params.ndraw

; initialize the output data structure
    result = {$
      chi2:               1E6,$ ; 1E6=not fitted
      mstar:    fltarr(ndraw),$
      age:      fltarr(ndraw),$
      sfrage:   fltarr(ndraw),$
      tau:      fltarr(ndraw),$
      zmetal:   fltarr(ndraw),$
      AV:       fltarr(ndraw),$
      mu:       fltarr(ndraw),$
      sfr:      fltarr(ndraw),$
      sfr100:   fltarr(ndraw),$
      b100:     fltarr(ndraw),$
      b1000:    fltarr(ndraw),$
      ewoii:    fltarr(ndraw),$
      ewoiiihb: fltarr(ndraw),$
      ewniiha:  fltarr(ndraw)}
    result = replicate(result,ngal)
    result.chi2 = isedfit.chi2

; optionally restore the spectra!

; build the posteriors
    good = where(result.chi2 lt 0.9E6,ngood)
    if ngood gt 0L then begin
       logscale_err = post[good].totalmass_err/post[good].totalmass/alog(10.0)
       logscale = alog10(post[good].totalmass)+randomn(seed,ndraw,ngal)*logscale_err

       result[good].mstar = alog10(modelgrid[post[good].draws].mstar)+logscale   ; [log Msun]
       result[good].sfr = modelgrid[post[good].draws].sfr+logscale               ; [Msun/yr]
       result[good].sfr100 = modelgrid[post[good].draws].sfr100+logscale         ; [Msun/yr]
       result[good].b100 = modelgrid[post[good].draws].b100
       result[good].b1000 = modelgrid[post[good].draws].b1000

       result[good].age = modelgrid[post[good].draws].age             ; [Gyr]
       result[good].sfrage = modelgrid[post[good].draws].sfrage       ; [Gyr]
       result[good].Zmetal = modelgrid[post[good].draws].Zmetal       ; [Zsun]
       result[good].tau = modelgrid[post[good].draws].tau             ; [Gyr]
       result[good].AV = modelgrid[post[good].draws].AV               ; [mag]
       result[good].mu = modelgrid[post[good].draws].mu
       result[good].ewoii = modelgrid[post[good].draws].ewoii         ; [Angstrom, rest]
       result[good].ewoiiihb = modelgrid[post[good].draws].ewoiiihb   ; [Angstrom, rest]
       result[good].ewniiha = modelgrid[post[good].draws].ewniiha     ; [Angstrom, rest]
    endif

; restore the photometric redshifts
    if tag_exist(post,'pofz') then begin
       result = struct_addtags(temporary(result),$
         struct_trimtags(post,select='pofz'))
    endif
    
return, result
end
