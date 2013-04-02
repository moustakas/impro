;+
; NAME:
;   ISEDFIT_RECONSTRUCT_POSTERIOR()
;
; PURPOSE:
;   Reconstruct the posterior distribution(s) of the output parameters.
;
; INPUTS:
;   paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   isedfit_dir - I/O path
;
; KEYWORD PARAMETERS:
;   maxold - see ISEDFIT
;   silent - suppress messages to STDOUT
;
; OUTPUTS:
;   model - data structure array containing the best-fitting spectrum
;     for each object (the structure will be different depending on
;     whether or not MAXOLD=1)  
;
; OPTIONAL OUTPUTS:
;   isedfit - ISEDFIT result structure
;
; COMMENTS:
;   The cosmological parameters are hard-wired to match
;   ISEDFIT_MODELS.  
;
;   Floating underflow messages are suppressed.
;
;   Note that the best-fitting models are interpolated onto a
;   common wavelength grid, currently hard-wired between 100 A and
;   5 microns, with a spacing of 2 A/pixel.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 27, NYU - largely excised from
;     ISEDFIT_QAPLOT and ISEDFIT_MEASURE
;
; Copyright (C) 2007, John Moustakas
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

function isedfit_reconstruct_posterior, isedfit_paramfile, post=post, params=params, $
  supergrid_paramfile=supergrid_paramfile, thissupergrid=thissupergrid, $
  isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, index=index, $
  outprefix=outprefix, age=age, sfrage=sfrage, tau=tau, Z=Z, av=av, nburst=nburst, $
  sfr0=sfr0, sfr100=sfr100, b100=b100, mgal=mgal, chunkindx=chunkindx, $
  modelindx=modelindx, indxage=ageindx, bigsfr0=bigsfr, bigmass=bigmass, bigsfrage=bigsfrage

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_reconstruct_posterior'
       return, -1
    endif

    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = './'
    if (n_elements(montegrids_dir) eq 0) then montegrids_dir = isedfit_dir+'montegrids/'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile)

; read the required SUPERGRID parameter file
    if n_elements(supergrid_paramfile) eq 0 then begin
       splog, 'SUPERGRID parameter file required'
       return, -1
    endif

    if n_elements(thissupergrid) eq 0 then begin
       splog, 'THISSUPERGRID must be specified!'
       return, -1
    endif
    super = read_supergrid_paramfile(supergrid_paramfile,supergrid=thissupergrid)

    fp = isedfit_filepaths(params,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=thissupergrid,isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
      outprefix=outprefix)

; restore the POSTERIOR output if not passes
    if (n_elements(post) eq 0L) then begin
       if (file_test(fp.isedfit_dir+fp.post_outfile+'.gz',/regular) eq 0L) then $
         message, 'POSTERIOR output not found!'
       splog, 'Reading '+fp.isedfit_dir+fp.post_outfile+'.gz'
       post = mrdfits(fp.isedfit_dir+fp.post_outfile+'.gz',1,/silent,rows=index)
    endif
    ngal = n_elements(post)

;   if (n_elements(isedfit) eq 0L) then begin
;      if (file_test(fp.isedfit_dir+fp.isedfit_outfile+'.gz',/regular) eq 0L) then $
;        message, 'ISEDFIT output not found!'
;      splog, 'Reading '+fp.isedfit_dir+fp.isedfit_outfile+'.gz'
;      isedfit = mrdfits(fp.isedfit_dir+fp.isedfit_outfile+'.gz',1,/silent,rows=index)
;   endif
    
; need all the chunks in memory!
    nchunk = n_elements(fp.sfhgrid_chunkfiles)
    for jj = 0, nchunk-1 do begin
       modelgrid1 = mrdfits(fp.modelspath+fp.isedfit_models_chunkfiles[jj]+'.gz',1,/silent)
       modelgrid1 = struct_trimtags(temporary(modelgrid1),except='modelmaggies')
       if (jj eq 0) then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].age)
    nallmodel = nmodel*nage
    ndraw = isedfit_ndraw() ; number of random draws

    bigmass = reform(modelgrid.mstar,nallmodel)
    bigage = reform(modelgrid.age,nallmodel)
    bigtau = reform(rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel),nallmodel)
    bigZ = reform(rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel),nallmodel)
    bigav = reform(rebin(reform(modelgrid.mu*modelgrid.av,1,nmodel),nage,nmodel),nallmodel)
    bignburst = reform(rebin(reform(modelgrid.nburst,1,nmodel),nage,nmodel),nallmodel)

; for restoring the posterior on the models    
    bigchunkindx = reform(rebin(reform(modelgrid.chunkindx,1,nmodel),nage,nmodel),nallmodel)
    bigmodelindx = reform(rebin(reform(modelgrid.modelindx,1,nmodel),nage,nmodel),nallmodel)
    bigageindx = reform(rebin(reform(lindgen(nage),nage,1),nage,nmodel),nallmodel)

    if arg_present(sfr0) or arg_present(sfr100) or arg_present(b100) or $
      arg_present(mgal) or arg_present(sfrage) then dosfr = 1 else dosfr = 0
    if dosfr then begin
       bigsfr = bigage*0D
       bigsfr100 = bigage*0D    ; average over the previous 100 Myr
       bigb100 = bigage*0D      ; birthrate parameter
       bigmgal = bigage*0D      ; galaxy mass ignoring mass loss 
       bigsfrage = bigage*0D    ; SFR-weighted age
       for imod = 0L, nmodel-1 do begin
          tindx = lindgen(nage)+imod*nage
          modelsfr = isedfit_reconstruct_sfh(modelgrid[imod],outage=bigage[tindx],$
            sfr100=modelsfr100,b100=modelb100,mgalaxy=modelmgal,sfrage=modelsfrage)

          bigsfr[tindx] = modelsfr       ; alog10(sfr)
          bigsfr100[tindx] = modelsfr100 ; alog10(sfr100) 
          bigb100[tindx] = modelb100
          bigmgal[tindx] = modelmgal
          bigsfrage[tindx] = modelsfrage
       endfor
; apply the scale factor
       b100 = fltarr(ndraw,ngal)
       sfr0 = fltarr(ndraw,ngal)
       sfr100 = fltarr(ndraw,ngal)
       sfrage = fltarr(ndraw,ngal)
;      mgal = fltarr(ndraw,ngal)
       for gg = 0L, ngal-1 do begin
          logscale_err = post[gg].scale_err/post[gg].scale/alog(10)
          logscale = alog10(post[gg].scale) + randomn(seed,ndraw)*logscale_err
          b100[*,gg] = alog10(bigb100[post[gg].draws])
          sfr0[*,gg] = alog10(bigsfr[post[gg].draws])+logscale
          sfr100[*,gg] = alog10(bigsfr100[post[gg].draws])+logscale
          sfrage[*,gg] = bigsfrage[post[gg].draws]
;         mgal[*,gg] = bigmgal[post[gg].draws]
       endfor    
    endif

; now get the remaining parameters    
    mass = fltarr(ndraw,ngal)
    for gg = 0L, ngal-1 do begin
       logscale_err = post[gg].scale_err/post[gg].scale/alog(10)
       logscale = alog10(post[gg].scale) + randomn(seed,ndraw)*logscale_err
       mass[*,gg] = alog10(bigmass[post[gg].draws])+logscale
    endfor

    if arg_present(age) then age = bigage[post.draws]
    if arg_present(tau) then tau = bigtau[post.draws]
    if arg_present(Z) then Z = bigZ[post.draws]
    if arg_present(av) then av = bigav[post.draws]
    if arg_present(chunkindx) then chunkindx = bigchunkindx[post.draws]
    if arg_present(modelindx) then modelindx = bigmodelindx[post.draws]
    if arg_present(ageindx) then ageindx = bigageindx[post.draws]
        
return, mass
end
