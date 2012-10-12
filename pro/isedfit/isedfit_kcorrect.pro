;+
; NAME:
;   ISEDFIT_KCORRECT
;
; PURPOSE:
;   Compute K-corrections from the best-fitting iSEDfit model. 
;
; INPUTS:
;   paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   isedpath - I/O path
;   index - zero-indexed list of objects to analyze (default is to
;     compute K-corrections for everything)
;   isedfit_sfhgrid_dir - see BUILD_ISEDFIT_SFHGRID
;   outprefix - optional output prefix string (see ISEDFIT) 
;   isedfit_sfhgrid_dir - see BUILD_ISEDFIT_SFHGRID
;   thissfhgrid - if SFHGRID in the parameter is a vector, then use
;     this optional input to pick one
;   galchunksize - divide the sample into GALCHUNKSIZE chunks to
;     minimize the amount of memory used (default 5000L)
;
;   out_filterlist - list of output filters in which outputs
;     (KCORRECT, ABSMAG, etc.) are desired (default SDSS ugriz)
;     [NFILT] 
;   band_shift - use band-shifted rest-frame bandpasses (default 0.0) 
;
; KEYWORD PARAMETERS:
;   vega - convert the output absolute magnitudes to Vega (default is
;     to keep them as AB magnitudes)
;
; OUTPUTS:
;   kcorrect - derived K-corrections [NFILT,NGAL]
;   absmag - absolute (rest-frame) magnitudes [NFILT,NGAL]
;   ivarabsmag - corresponding inverse variance [NFILT,NGAL]
;   synth_absmag - like ABSMAG, but as synthesized from the
;     best-fitting model [NFILT,NGAL]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is a glorified wrapper on IM_SIMPLE_KCORRECT. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Aug 30, UCSD
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

pro isedfit_kcorrect, paramfile, isedfit, params=params, isedpath=isedpath, index=index, $
  outprefix=outprefix, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, thissfhgrid=thissfhgrid, $
  galchunksize=galchunksize, out_filterlst=out_filterlist, band_shift=band_shift, $
  kcorrect=kcorrect, absmag=absmag, ivarabsmag=ivarabsmag, synth_absmag=synth_absmag, $
  vega=vega
    
    if (n_elements(paramfile) eq 0L) and (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_kcorrect'
       return
    endif

; read the parameter file and parse to get the relevant path and
; filenames  
    if (n_elements(isedpath) eq 0) then isedpath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)

; make sure SFHGRID is a scalar
    nsfhgrid = n_elements(params.sfhgrid)
    if (nsfhgrid gt 1) then begin
       if (n_elements(thissfhgrid) eq 0) then begin
          splog, 'SFHGRID in the parameter file is a vector!'
          splog, 'Use the THISSFHGRID optional input to choose one'
          return
       endif else begin
          params = struct_trimtags(params,except='sfhgrid')
          params = struct_addtags(params,{sfhgrid: thissfhgrid})
       endelse
    endif

    fp = isedfit_filepaths(params,outprefix=outprefix,isedpath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    kcorrfile = isedpath+strtrim(fp.kcorr_outfile,2)
    if file_test(kcorrfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+kcorrfile+' exists; use /CLOBBER'
       return
    endif

; restore the filters and the results
    filterlist = strtrim(params.filterlist,2)

    nfilt = n_elements(out_filterlist)
    if (nfilt eq 0) then begin
       splog, 'No output filters defined! Adopting ugriz'
       out_filterlist = sdss_filterlist()
    endif
    
    junk = isedfit_restore(paramfile,isedfit,params=params,$
      isedpath=isedpath,index=index,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,$
      outprefix=outprefix,silent=silent,/nomodels)
    ngal = n_elements(isedfit)

; compute K-corrections; split the problem into chunks because the
; arrays can be memory intensive for large samples
    if (n_elements(galchunksize) eq 0) then galchunksize = 5000L
    ngalchunk = ceil(ngal/float(galchunksize))
    sortindx = sort(isedfit.chunkindx) ; sort for speed

    kcorrect = fltarr(nfilt,ngal)
    absmag = fltarr(nfilt,ngal)
    ivarabsmag = fltarr(nfilt,ngal)
    synth_absmag = fltarr(nfilt,ngal)
    
    t0 = systime(1)
    for gchunk = 0L, ngalchunk-1L do begin
       t1 = systime(1)
       splog, format='("Computing K-corrections for galaxy chunk ",I0,"/",I0)', $
         gchunk+1, ngalchunk

       g1 = gchunk*galchunksize
       g2 = ((gchunk*galchunksize+galchunksize)<ngal)-1L
       gnthese = g2-g1+1L & gthese = lindgen(gnthese)+g1

       good = where((isedfit[sortindx[gthese]].chi2 gt 0.0) and $
         (isedfit[sortindx[gthese]].chi2 lt 1E6),ngood)
       if (ngood ne 0L) then begin
          model = isedfit_restore(paramfile,params=params,isedpath=isedpath,$
            index=sortindx[gthese[good]],/flambda,outprefix=outprefix,$
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,silent=silent)
          oneplusz = rebin(reform(1.0+isedfit[sortindx[gthese[good]]].zobj,1,ngood),$
            n_elements(model[0].wave),ngood)
          restwave = model.wave/oneplusz
          restflux = model.flux*oneplusz

          chunk_kcorr = im_simple_kcorrect(isedfit[sortindx[gthese[good]]].zobj,$
            isedfit[sortindx[gthese[good]]].maggies,isedfit[sortindx[gthese[good]]].ivarmaggies,$
            filterlist,out_filterlist,restwave,restflux,band_shift=band_shift,$
            absmag=chunk_absmag,ivarabsmag=chunk_ivarabsmag,synth_absmag=chunk_synth_absmag,$
            vega=vega,/silent)

          kcorrect[*,sortindx[gthese[good]]] = chunk_kcorr
          absmag[*,sortindx[gthese[good]]] = chunk_absmag
          ivarabsmag[*,sortindx[gthese[good]]] = chunk_ivarabsmag
          synth_absmag[*,sortindx[gthese[good]]] = chunk_synth_absmag
       endif 
       if (gchunk eq 0) then splog, format='("Time for first '+$
         'Chunk = ",G0," minutes          ")', (systime(1)-t1)/60.0
    endfor
    splog, format='("Time for all Chunks = ",G0," minutes          ")', $
      (systime(1)-t0)/60.0
    
return
end
