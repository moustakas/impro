;+
; NAME:
;   ISEDFIT_KCORRECT
;
; PURPOSE:
;   Measure various quantities from the best-fitting
;   star-formation histories derived by ISEDFIT.
;
; INPUTS/OUTPUTS:
;   result      - see ISEDFIT
;   result_info - see ISEDFIT
;
; OPTIONAL INPUTS:
;   datapath       - should match ISEDFIT [default CWD()] 
;   isedfitprefix  - should match ISEDFIT (default 'isedfit')
;   restfilterfile - see COMMENTS, below
;
; KEYWORD PARAMETERS:
;   maxold - see ISEDFIT
;   write  - write out
;
; OUTPUTS:
;   measure - output data structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   RESTFILTERFILE can be used to specify the rest-frame
;   magnitudes of interest other than the defaults.  The file
;   should contain three columns containing the filter name (which
;   must be in the KCORRECT database), the bandpass name (which
;   will be the output structure tag name), and a Boolean flag
;   indicating whether the output magnitude should be in Vega (1)
;   or AB (0).  For example, the file might look like this:
;
;     bessell_B.par B       1
;     bessell_V.par V       1
;     sdss_g0.par   sdss_g  0
;     twomass_J.par J       1
;
; TODO: 
;   (1) Compute stellar mass-to-light ratios.
;   (2) Add more error checking when reading RESTFILTERFILE. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 March 07, U of A - written
;   jm06aug02uofa - major re-write
;   jm07mar21nyu  - further developments; speeded up 
;   jm07jun28nyu  - updated 
;
; Copyright (C) 2006-2007, John Moustakas
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

pro isedfit_kcorrect, paramfile, isedfit, params=params, iopath=iopath, $
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
    if (n_elements(iopath) eq 0) then iopath = './'
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

    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    kcorrfile = iopath+strtrim(fp.kcorr_outfile,2)
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
      iopath=iopath,index=index,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,$
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
          model = isedfit_restore(paramfile,params=params,iopath=iopath,$
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
