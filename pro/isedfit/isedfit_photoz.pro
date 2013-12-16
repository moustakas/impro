;+
; NAME:
;   ISEDFIT_PHOTOZ
;
; PURPOSE:
;   Photometric-redshift wrapper on iSEDfit to be used when the
;   redshift is unknown (or only approximately known).
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;   maggies - input galaxy photometry [NFILT,NGAL]
;   ivarmaggies - inverse variance array for MAGGIES [NFILT,NGAL]  
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the results should be
;     written (default PWD=present working directory)  
;
;   outprefix - optionally write out files with a different prefix
;     from that specified in ISEDFIT_PARAMFILE (or PARAMS) (very
;     useful for fitting various subsets of the input sample with
;     different assumptions)
;   index - use this optional input to fit a zero-indexed subset of
;     the full sample (default is to fit everything)
;
;   ra,dec - galaxy right ascension and declination which will be
;     copied to the output structure for the user's convenience
;     [NGAL] (decimal degrees)
;
; KEYWORD PARAMETERS:
;   silent - suppress messages to STDOUT
;   nowrite - do not write out any of the output files (generally not
;     recommended but can be useful in certain situations) 
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;   extra - additional keywords to pass to ISEDFIT
;
; OUTPUTS:
;   Binary FITS tables containing the fitting results are written
;   to ISEDFIT_DIR. 
;
; OPTIONAL OUTPUTS:
;   isedfit_photoz_results - output data structure containing all the
;     results; see the iSEDfit documentation for a detailed breakdown
;     and explanation of all the outputs
;   isedfit_photoz_post - output data structure containing the random draws
;     from the posterior distribution function, which can be used to
;     rebuild the posterior distributions of any of the output
;     parameters (using ISEDFIT_RECONSTRUCT_POSTERIOR) 
;   isedfit_photoz_outfile - output file name of the iSEDfit photoz
;     results 
;
; COMMENTS:
;   Better documentation of the output data structures would be
;   helpful here.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Dec 04, Siena
;
; Copyright (C) 2013, John Moustakas
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

pro isedfit_photoz, isedfit_paramfile, maggies, ivarmaggies, params=params, $
  thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, outprefix=outprefix, $
  index=index, ra=ra, dec=dec, isedfit_photoz_results=isedfit_photoz_results, $
  isedfit_photoz_post=isedfit_photoz_post, isedfit_photoz_outfile=isedfit_photoz_outfile, $
  silent=silent, nowrite=nowrite, clobber=clobber, _extra=extra

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_photoz'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()

; error checking on the input photometry    
    ndim = size(maggies,/n_dim)
    dims = size(maggies,/dim)
    if (ndim eq 1) then ngal = 1 else ngal = dims[1]  ; number of galaxies
    nfilt = dims[0] ; number of filters

    nmaggies = n_elements(maggies)
    nivarmaggies = n_elements(ivarmaggies)
    
    if (nmaggies eq 0L) or (nivarmaggies eq 0L) then begin
       doc_library, 'isedfit_photoz'
       return
    endif

    ndim = size(maggies,/n_dimension)
    if (ndim ne 2L) then begin ; reform into a 2D array
       maggies = reform(maggies,n_elements(maggies),1)
       ivarmaggies = reform(ivarmaggies,n_elements(maggies),1)
    endif

    if (n_elements(maggies) ne n_elements(ivarmaggies)) then begin
       splog, 'Dimensions of MAGGIES and IVARMAGGIES do not match'
       return
    endif
    if (total(finite(maggies) eq 0B) ne 0.0) or $
      (total(finite(ivarmaggies) eq 0B) ne 0.0) then begin
       splog, 'MAGGIES and IVARMAGGIES cannot have infinite values!'
       return
    endif

; check for RA,DEC    
    if n_elements(ra) ne 0L then begin
       if size(ra,/type) ne 5 then splog, 'Warning: RA should be type double!'
       if n_elements(ra) ne ngal then begin
          splog, 'Dimensions of RA must match the input number of galaxies!'
          return
       endif
    endif else ra = dblarr(ngal)-1D

    if n_elements(dec) ne 0L then begin
       if size(dec,/type) ne 5 then splog, 'Warning: DEC should be type double!'
       if n_elements(dec) ne ngal then begin
          splog, 'Dimensions of DEC must match the input number of galaxies!'
          return
       endif
    endif else dec = dblarr(ngal)-1D
    
; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          isedfit_photoz, isedfit_paramfile1, maggies, ivarmaggies, params=params[ii], $
            isedfit_dir=isedfit_dir, outprefix=outprefix, index=index, $
            ra=ra, dec=dec, isedfit_photoz_results=isedfit_photoz_results, $
            isedfit_photoz_post=isedfit_photoz_post, isedfit_photoz_outfile=isedfit_photoz_outfile, $
            silent=silent, nowrite=nowrite, clobber=clobber, _extra=extra
       endfor 
       return
    endif

    fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,outprefix=outprefix)

; hack!!    use the regular _outfile
    isedfit_photoz_outfile = fp.isedfit_dir+fp.isedfit_outfile
;   isedfit_photoz_outfile = fp.isedfit_dir+fp.isedfit_photoz_outfile
    if file_test(isedfit_photoz_outfile+'.gz') and (keyword_set(clobber) eq 0) and $
      (keyword_set(nowrite) eq 0) then begin
       splog, 'Output file '+isedfit_photoz_outfile+' exists; use /CLOBBER'
       return
    endif

;; fit the requested subset of objects and return
;    if (n_elements(index) ne 0L) then begin
;       isedfit, isedfit_paramfile1, maggies[*,index], ivarmaggies[*,index], $
;         z[index], params=params, isedfit_dir=isedfit_dir, outprefix=outprefix, $
;         ra=ra[index], dec=dec[index], isedfit_results=isedfit_results1, $
;         isedfit_post=isedfit_post1, isedfit_outfile=isedfit_outfile, allages=allages, $
;         maxold=maxold, silent=silent, /nowrite, clobber=clobber
;
;       isedfit_results = init_isedfit(ngal,nfilt,params=params,$
;         ra=ra,dec=dec,isedfit_post=isedfit_post)
;       isedfit_results[index] = isedfit_results1
;       isedfit_post[index] = isedfit_post1
;       if (keyword_set(nowrite) eq 0) then begin
;          im_mwrfits, isedfit_results, isedfit_outfile, /clobber, silent=silent
;          im_mwrfits, isedfit_post, fp.isedfit_dir+fp.post_outfile, /clobber, silent=silent
;       endif
;       return
;    endif
    
; trick iSEDfit into thinking we have a bigger sample than we actually
; do; the array juggling here might barf if the dimensions get too
; large but punt on dealing with this for now
    redshift = params.redshift
    nz = n_elements(redshift)
    nfilt = n_elements(params.filterlist)

    maggies1 = reform(rebin(reform(maggies,nfilt,ngal,1),nfilt,ngal,nz),nfilt,ngal*nz)
    ivarmaggies1 = reform(rebin(reform(ivarmaggies,nfilt,ngal,1),nfilt,ngal,nz),nfilt,ngal*nz)
    redshift1 = reform(rebin(reform(redshift,1,nz),ngal,nz),ngal*nz)
    
    isedfit, isedfit_paramfile, maggies1, ivarmaggies1, redshift1, $
      params=params, isedfit_dir=isedfit_dir, isedfit_results=ised, $
      isedfit_post=isedpost, /nowrite

; need some new code here to deal with the posteriors; ISEDFIT_POST
; probably won't work

; build the output data structure
    ised = reform(ised,ngal,nz)
    isedfit_photoz_results = struct_trimtags(im_empty_structure(ised,ncopies=ngal),$
      select=['isedfit_id','ra','dec','maggies','ivarmaggies','bestmaggies',$
      'chi2','mstar','age','sfrage','tau','zmetal','av','mu','oiiihb',$
      'nlyc','sfr','sfr100','b100','b1000','ewoii','ewoiiihb','ewniiha'])
    isedfit_photoz_results = struct_addtags(isedfit_photoz_results,replicate($
      {redshift: redshift, pofz: fltarr(nz), chi2_z: 0.0, z: 0.0, z_50: 0.0, $
      z_avg: 0.0, z_err: 0.0},ngal))
    for igal = 0, ngal-1 do begin
       pofz = exp(-0.5*(ised[igal,*].chi2-min(ised[igal,*].chi2)))
       ml = max(pofz,mindx)
       isedfit_photoz_results[igal] = im_struct_assign(ised[igal,mindx],isedfit_photoz_results[igal],/nozer)
       isedfit_photoz_results[igal].isedfit_id = igal

       isedfit_photoz_results[igal].pofz = pofz
       isedfit_photoz_results[igal].chi2_z = ised[igal,mindx].chi2
       isedfit_photoz_results[igal].z = redshift[mindx]
;      isedfit_photoz_results[igal].z_50 = 
;      isedfit_photoz_results[igal].z_avg = 
;      isedfit_photoz_results[igal].z_err = 
    endfor
    
; write out the final structure and the full posterior distributions
    if keyword_set(nowrite) eq 0 then begin
       im_mwrfits, isedfit_photoz_results, isedfit_photoz_outfile, silent=silent, /clobber
;      im_mwrfits, isedfit_photoz_post, fp.isedfit_dir+fp.post_photoz_outfile, silent=silent, /clobber
    endif

return
end
