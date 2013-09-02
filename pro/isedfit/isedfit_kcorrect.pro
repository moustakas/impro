;+
; NAME:
;   ISEDFIT_KCORRECT
;
; PURPOSE:
;   Compute K-corrections from the best-fitting iSEDfit model.
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where all the ISEDFIT output
;     files can be found; should match the directory passed to ISEDFIT
;     (default PWD=present working directory)  
;   montegrids_dir - full directory path where the Monte Carlo grids
;     written by ISEDFIT_MONTEGRIDS can be found (default 'montegrids'
;     subdirectory of the PWD=present working directory)
;   index - zero-indexed list of objects to analyze (default is to
;     compute K-corrections for everything)
;   outprefix - optional output prefix string (see ISEDFIT) 
;
;   absmag_filterlist - list of output filters in which K-corrected
;     absolute magnitudes are desired (default SDSS ugriz) [NFILT] 
;   band_shift - use band-shifted rest-frame bandpasses (default 0.0) 
;
; KEYWORD PARAMETERS:
;   vega - convert the output absolute magnitudes to Vega (default is
;     to keep them as AB magnitudes)
;   nowrite - do not write out any of the output files (generally not
;     recommended but can be useful in certain situations) 
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;
; OUTPUTS:
;   kcorrect_results - [NGAL] data structure with the following tags
;     in addition to some handy info culled from the ISEDFIT data
;     structure: 
;      .kcorrect - derived K-corrections [NFILT]
;      .absmag - absolute (rest-frame) magnitudes [NFILT]
;      .ivarabsmag - corresponding inverse variance [NFILT]
;      .synth_absmag - like ABSMAG, but as synthesized from the
;        best-fitting model [NFILT] 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is a glorified wrapper on IM_SIMPLE_KCORRECT. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 30, UCSD
;   jm13aug09siena - updated to conform to the latest data model;
;     documentation updated 
;
; Copyright (C) 2011, 2013, John Moustakas
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

pro isedfit_kcorrect, isedfit_paramfile, params=params, thissfhgrid=thissfhgrid, $
  isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, outprefix=outprefix, $
  index=index, absmag_filterlist=absmag_filterlist, band_shift=band_shift, $
  kcorrect_results=kcorrect_results, vega=vega, nowrite=nowrite, clobber=clobber

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_kcorrect'
       return
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
          isedfit_kcorrect, params=params[ii], isedfit_dir=isedfit_dir,$
            montegrids_dir=montegrids_dir, outprefix=outprefix, index=index, $
            absmag_filterlist=absmag_filterlist, band_shift=band_shift, $
            kcorrect_results=kcorrect_results1, vega=vega, nowrite=nowrite, $
            clobber=clobber
          if ii eq 0 then kcorrect_results = kcorrect_results1 else $
            kcorrect_results = [[temporary(kcorrect_results)],[kcorrect_results1]]
       endfor 
       return
    endif 
       
    fp = isedfit_filepaths(params,outprefix=outprefix,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,band_shift=band_shift)
    kcorrfile = fp.isedfit_dir+fp.kcorr_outfile
    if file_test(kcorrfile) and keyword_set(clobber) eq 0 then begin
       splog, 'Output file '+kcorrfile+' exists; use /CLOBBER'
       return
    endif

; input/output filters
    filterlist = strtrim(params.filterlist,2)
    nfilt = n_elements(absmag_filterlist)
    if (nfilt eq 0) then begin
       splog, 'No output filters defined! Adopting SDSS/ugriz'
       absmag_filterlist = sdss_filterlist()
       nfilt = n_elements(absmag_filterlist)
    endif

; restore the ISEDFIT results without the best-fitting models to get
; the number of galaxies and other goodies
    isedfit = read_isedfit(params=params,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,index=index,outprefix=outprefix,$
      silent=silent)
    ngal = n_elements(isedfit)

; initialize the output data structure    
    kcorrect_results = struct_trimtags(isedfit,select=['isedfit_id',$
      'z','maggies','ivarmaggies','chi2'])
    kcorrect_results = struct_addtags(temporary(kcorrect_results),replicate({$
      kcorrect:          fltarr(nfilt), $
      absmag:            fltarr(nfilt),$
      ivarabsmag:        fltarr(nfilt),$
      synth_absmag:      fltarr(nfilt),$
      absmag_filterlist: absmag_filterlist},ngal))
    
; compute K-corrections; split the problem into chunks because the
; arrays can be memory intensive for large samples
    ngalchunk = ceil(ngal/float(params.galchunksize))
    sortindx = sort(isedfit.chunkindx) ; sort for speed

    t0 = systime(1)
    for gchunk = 0L, ngalchunk-1 do begin
       t1 = systime(1)
       splog, format='("Computing K-corrections for galaxy chunk ",I0,"/",I0)', $
         gchunk+1, ngalchunk

       g1 = gchunk*params.galchunksize
       g2 = ((gchunk*params.galchunksize+params.galchunksize)<ngal)-1L
       gnthese = g2-g1+1L 
       gthese = lindgen(gnthese)+g1

       good = where((isedfit[sortindx[gthese]].chi2 gt 0.0) and $
         (isedfit[sortindx[gthese]].chi2 lt 0.9E6),ngood)
       if (ngood ne 0L) then begin
          isedfit1 = read_isedfit(params=params,isedfit_dir=isedfit_dir,$
            montegrids_dir=montegrids_dir,index=index[sortindx[gthese[good]]],$
            outprefix=outprefix,/flambda,/silent,/getmodels)
          oneplusz = rebin(reform(1.0+isedfit[sortindx[gthese[good]]].z,1,ngood),$
            n_elements(isedfit1[0].wave),ngood)
          restwave = isedfit1.wave/oneplusz
          restflux = isedfit1.flux*oneplusz

          chunk_kcorr = im_simple_kcorrect(isedfit[sortindx[gthese[good]]].z,$
            isedfit[sortindx[gthese[good]]].maggies,isedfit[sortindx[gthese[good]]].ivarmaggies,$
            filterlist,absmag_filterlist,restwave,restflux,band_shift=band_shift,$
            absmag=chunk_absmag,ivarabsmag=chunk_ivarabsmag,synth_absmag=chunk_synth_absmag,$
            chi2=chi2,vega=vega,/silent)

          kcorrect_results[sortindx[gthese[good]]].kcorrect = chunk_kcorr
          kcorrect_results[sortindx[gthese[good]]].absmag = chunk_absmag
          kcorrect_results[sortindx[gthese[good]]].ivarabsmag = chunk_ivarabsmag
          kcorrect_results[sortindx[gthese[good]]].synth_absmag = chunk_synth_absmag
       endif 
       if (gchunk eq 0) then splog, format='("Time for first '+$
         'Chunk = ",G0," minutes          ")', (systime(1)-t1)/60.0
    endfor
    splog, format='("Time for all Chunks = ",G0," minutes          ")', $
      (systime(1)-t0)/60.0

; write out
    if keyword_set(nowrite) eq 0 then im_mwrfits, kcorrect_results, $
      kcorrfile, silent=silent, /clobber
    
return
end
