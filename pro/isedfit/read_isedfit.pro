;+
; NAME:
;   READ_ISEDFIT()
;
; PURPOSE:
;   Restore the ISEDFIT fitting results and, optionally, the posterior
;   parameter distributions, and the best-fitting models. 
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where all the ISEDFIT output
;     files can be found; should match the directory passed to ISEDFIT
;     (default PWD=present working directory)  
;   montegrids_dir - full directory path where the Monte Carlo grids
;     written by ISEDFIT_MONTEGRIDS can be found (default 'montegrids'
;     subdirectory of the PWD=present working directory)
;   index - use this optional input to restore a zero-indexed subset
;     of the full sample (default is to restore everything, although
;     see COMMENTS)
;   outprefix - optional output prefix string (see ISEDFIT) 
;   in_isedfit - use this input data structure in lieu of the ISEDFIT
;     output structure read based on the contents of
;     ISEDFIT_PARAMFILE; this input facilitates some of the
;     fancy plotting I've done for papers  
;
; KEYWORD PARAMETERS:
;   flambda - convert to F(lambda) (erg/s/cm^2/A) (default is AB mag)
;   fnu - convert to F(nu) (erg/s/cm^2/Hz) (default is AB mag) 
;   getmodels - restore the ISEDFIT output structure *and* the model
;     spectra 
;   noigm - do not convolve with the IGM, over-riding the content of
;     ISEDFIT_PARAMFILE; this can be useful for testing but should not
;     in general be used
;   silent - suppress messages to STDOUT
;
; OUTPUTS:
;   model - output data structure array containing all the ISEDFIT
;     outputs plus the following tags if /GETMODELS [NGAL]
;     .WAVE - observed-frame wavelength array (Angstrom)
;     .FLUX - observed-frame flux array (AB mag)
;
; OPTIONAL OUTPUTS:
;   isedfit_post - iSEDfit result structure (see ISEDFIT) [NGAL] 
;
; COMMENTS:
;   This routine should be not be used to restore too many galaxy
;   spectra (if /GETMODELS) or posteriors (if ISEDFIT_POST is
;   requested), otherwise memory problems may occur; use INDEX. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 27, NYU - largely excised from
;     ISEDFIT_SED_QAPLOT and ISEDFIT_MEASURE
;   jm11oct01ucsd - documentation cleaned up and updated
;   jm13aug09siena - updated to conform to the latest data model;
;     documentation updated 
;
; Copyright (C) 2007, 2011, 2013, John Moustakas
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

function read_isedfit, isedfit_paramfile, params=params, thissfhgrid=thissfhgrid, $
  isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, in_isedfit=in_isedfit, $
  outprefix=outprefix, index=index, isedfit_post=isedfit_post, flambda=flambda, $
  fnu=fnu, getmodels=getmodels, noigm=noigm, silent=silent

    common isedfit_igmgrid, igmgrid
    
    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'read_isedfit'
       return, -1
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()
    if (n_elements(montegrids_dir) eq 0) then montegrids_dir = get_pwd()+'montegrids/'

; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
; there's no clean way (that I can think of!) to prevent
; isedfit_post from being read even if the user does not want it 
          if arg_present(isedfit_post) then begin
             result1 = read_isedfit(params=params[ii],isedfit_dir=isedfit_dir,$
               montegrids_dir=montegrids_dir,in_isedfit=in_isedfit,$
               outprefix=outprefix,index=index,isedfit_post=isedfit_post1,$
               flambda=flambda,fnu=fnu,getmodels=getmodels,noigm=noigm,silent=silent)
          endif else begin
             result1 = read_isedfit(params=params[ii],isedfit_dir=isedfit_dir,$
               montegrids_dir=montegrids_dir,in_isedfit=in_isedfit,$
               outprefix=outprefix,index=index,$;isedfit_post=isedfit_post1,$
               flambda=flambda,fnu=fnu,getmodels=getmodels,noigm=noigm,silent=silent)
          endelse
          if ii eq 0 then result = result1 else result = [[result],[result1]]
          if arg_present(isedfit_post) then begin
             if ii eq 0 then isedfit_post = temporary(isedfit_post1) else $
               isedfit_post = [[temporary(isedfit_post)],[temporary(isedfit_post1)]]
          endif
       endfor
       return, result
    endif 
       
; restore the ISEDFIT output; optionally take ISEDFIT as an input
; structure
    fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,outprefix=outprefix)
    if (n_elements(in_isedfit) eq 0L) then begin
       isedfile = fp.isedfit_dir+fp.isedfit_outfile+'.gz'
       if file_test(isedfile) eq 0 then begin
          splog, 'iSEDfit output file '+isedfile+' not found!'
          return, -1
       endif
       if keyword_set(silent) eq 0 then splog, 'Reading '+isedfile
       isedfit1 = mrdfits(isedfile,1,/silent)
; only restore a subset of the objects    
       if (n_elements(index) eq 0L) then begin
          index = lindgen(n_elements(isedfit1)) ; need this for isedfit_kcorrect
          isedfit = temporary(isedfit1) 
       endif else isedfit = temporary(isedfit1[index])
    endif else isedfit = in_isedfit
    ngal = n_elements(isedfit)

; restore the posteriors if requested
    if arg_present(isedfit_post) then begin
       isedfit_post = isedfit_reconstruct_posterior(params=params,$
         isedfit_dir=isedfit_dir,index=index,outprefix=outprefix)
    endif

; restore the models    
    if keyword_set(getmodels) eq 0 then begin
       return, isedfit ; done!
    endif else begin
       light = 2.99792458D18             ; speed of light [A/s]
       pc10 = 3.085678D19                ; fiducial distance [10 pc in cm]
       dlum = pc10*10D^(lf_distmod(isedfit.z,omega0=params.omega0,$ ; [cm]
         omegal0=params.omegal)/5D)/params.h100 
;      dlum = dluminosity(isedfit.z,/cm) ; luminosity distance [cm]
       
; read the models; group by chunks; remove objects that were not
; fitted
       allchunks = isedfit.chunkindx
       chunks = allchunks[uniq(allchunks,sort(allchunks))]
       nchunk = n_elements(chunks)
       
; initialize the model structure (need to get the number of pixels)
       montefile = strtrim(fp.montegrids_chunkfiles[0],2) ; check the first file
       if file_test(montefile+'*') eq 0 then begin
          splog, 'First ISEDFIT_MONTEGRIDS ChunkFile '+montefile+' not found!'
          return, -1
       endif

       junk = gz_mrdfits(fp.montegrids_chunkfiles[0],1,row=0,/silent)
       npix = n_elements(junk.wave)

       result = struct_addtags(temporary(isedfit),replicate($
         {wave: fltarr(npix), flux: fltarr(npix)},ngal))
       
       for ichunk = 0L, nchunk-1L do begin
          these = where(chunks[ichunk] eq allchunks,nthese)
          if (nthese ne 0L) and (chunks[ichunk] ge 0) then begin
             chunkfile = strtrim(fp.montegrids_chunkfiles[chunks[ichunk]],2)
             if (not keyword_set(silent)) then splog, 'Reading '+chunkfile
             grid = gz_mrdfits(chunkfile,1,/silent)
             for ii = 0L, nthese-1L do begin
                if (result[these[ii]].chi2 lt 1E6) then begin
                   i1 = result[these[ii]].modelindx ; mod 100
                   result[these[ii]].wave = grid[i1].wave
                   result[these[ii]].flux = grid[i1].flux
                endif
             endfor
          endif  
       endfor

; read the IGM absorption table
       if n_elements(igmgrid) eq 0 then begin
          if params.igm and (keyword_set(noigm) eq 0) then begin
             splog, 'Reading IGM attenuation lookup table'
             igmgrid = mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)
          endif
       endif

; loop through every object and construct the best-fit model
       for igal = 0L, ngal-1L do begin
          if (result[igal].chi2 lt 1E6) then begin
             z = result[igal].z
             zwave = result[igal].wave*(1+z)
             zflux_flam = result[igal].totalmass*result[igal].flux*$ ; [erg/s/cm2/A]
               (pc10/dlum[igal])^2.0/(1.0+z)
             if params.igm and (keyword_set(noigm) eq 0) then begin
                windx = findex(igmgrid.wave,zwave)
                zindx = findex(igmgrid.zgrid,z)
                igm = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
                zflux_flam = zflux_flam*igm
             endif
             zflux_fnu = zflux_flam*zwave^2/light ; [erg/s/cm2/Hz]
             zflux_ab = -2.5*alog10(zflux_fnu>1D-50)-48.6
             result[igal].wave = zwave
             result[igal].flux = zflux_ab ; default
             if keyword_set(fnu) then result[igal].flux = zflux_fnu
             if keyword_set(flambda) then result[igal].flux = zflux_flam
          endif
       endfor
    endelse

return, result
end
