;+
; NAME:
;   ISEDFIT_MEASURE
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

function init_isedfit_measure, ngal, nfilt
; internal support routine - initialize the output data structure

    measure = {$
      isedfit_id:                       -1L, $
      zobj:                          -999.0, $
      maggies:          fltarr(nfilt)-999.0, $
      ivarmaggies:      fltarr(nfilt)-999.0, $
      bestmaggies:      fltarr(nfilt)-999.0, $
      chi2:                          -999.0, $
      galex_absmag:         fltarr(2)-999.0, $
      galex_absmag_ivar:    fltarr(2)-999.0, $
      galex_kcorrect:       fltarr(2)-999.0, $
      ugriz_absmag:         fltarr(5)-999.0, $
      ugriz_absmag_ivar:    fltarr(5)-999.0, $
      ugriz_kcorrect:       fltarr(5)-999.0, $
      ubvrijhk_absmag:      fltarr(8)-999.0, $
      ubvrijhk_absmag_ivar: fltarr(8)-999.0, $
      ubvrijhk_kcorrect:    fltarr(8)-999.0, $
;     select_absmag:                 -999.0, $
      cflux:                fltarr(5)-999.0, $ ; see IM_SIMPLE_KCORRECT
      uvflux:               fltarr(2)-999.0, $ ; see IM_SIMPLE_KCORRECT
      d4000_narrow:                  -999.0, $
      lick_hd_a:                     -999.0}
    measure = replicate(measure,ngal)

return, measure
end

pro isedfit_measure, paramfile, measure, isedfit, params=params, $
  iopath=iopath, outprefix=outprefix, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
  clobber=clobber, nowrite=nowrite, abmag=abmag, _extra=extra
; ABMAG - convert the UBVRIJHK rest-frame photometry to AB
    
    light = 2.99792458D18       ; speed of light [A/s]

    if (n_elements(paramfile) eq 0L) and (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_measure'
       return
    endif

; read the parameter file and parse to get the relevant path and
; filenames  
    if (n_elements(iopath) eq 0) then iopath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(iopath+paramfile)

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
             isedfit_measure, params=newparams2, iopath=iopath, $
               outprefix=outprefix, clobber=clobber, $
               nowrite=nowrite, abmag=abmag, _extra=extra
          endfor
       endfor 
       return 
    endif 

    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    measurefile = iopath+strtrim(fp.measure_outfile,2)
    if file_test(measurefile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+measurefile+' exists; use /CLOBBER'
       return
    endif

; restore the filters and the results
    filterlist = strtrim(params.filterlist,2)
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filterlist)

    junk = isedfit_restore(paramfile,isedfit,params=params,$
      iopath=iopath,outprefix=outprefix,silent=silent,/nomodels)
    ngal = n_elements(isedfit)

;   oneplusz = rebin(reform(1.0+measure[gthese[good]].zobj,1,ngal),$
;     n_elements(model[0].wave),ngal)
;   restwave = model.wave/oneplusz
;   restflux = model.flux*oneplusz

; initialize the output data structure 
    measure = init_isedfit_measure(ngal,nfilt)
    struct_assign, isedfit, measure, /nozero

;; compute spectral properties - too slow for large samples!!
;    good = where(isedfit.chi2 lt 1E6,ngood)
;    if (ngood ne 0L) then begin
;       for igal = 0L, ngood-1L do begin
;          spec = spectral_indices(restwave[*,good[igal]],$
;            restflux[*,good[igal]],/silent)
;          measure[good[igal]].d4000_narrow = spec.d4000_narrow[0]
;          measure[good[igal]].lick_hd_a = spec.lick_hd_a[0]
;       endfor
;    endif

; compute K-corrections; split the problem into chunks because the
; arrays can be memory intensive for large samples
    splog, 'Computing K-corrections'

    ubvrijhk_filterlist = [bessell_filterlist(),twomass_filterlist()]
    ugriz_filterlist = sdss_filterlist()

    sortindx = sort(isedfit.chunkindx) ; sort for speed

    galchunksize = 5000L
    ngalchunk = ceil(ngal/float(galchunksize))
    t0 = systime(1)
    for gchunk = 0L, ngalchunk-1L do begin
       t1 = systime(1)
       splog, format='("Chunk ",I0,"/",I0)', gchunk+1L, ngalchunk

       g1 = gchunk*galchunksize
       g2 = ((gchunk*galchunksize+galchunksize)<ngal)-1L
       gnthese = g2-g1+1L & gthese = lindgen(gnthese)+g1

       good = where((isedfit[sortindx[gthese]].chi2 gt 0.0) and $
         (isedfit[sortindx[gthese]].chi2 lt 1E6),ngood)
       if (ngood ne 0) then begin
          model = isedfit_restore(paramfile,params=params,iopath=iopath,$
            outprefix=outprefix,index=sortindx[gthese[good]],/flambda,silent=silent)
          oneplusz = rebin(reform(1.0+measure[sortindx[gthese[good]]].zobj,1,ngood),$
            n_elements(model[0].wave),ngood)
          restwave = model.wave/oneplusz
          restflux = model.flux*oneplusz

; UBVRIJHKs, band_shift=0.0, Vega       
          kcorr = im_simple_kcorrect(measure[sortindx[gthese[good]]].zobj,$
            measure[sortindx[gthese[good]]].maggies,measure[sortindx[gthese[good]]].ivarmaggies,$
            filterlist,ubvrijhk_filterlist,restwave,restflux,absmag=absmag,$
            ivarabsmag=absmag_ivar,clineflux=cflux,uvflux=uvflux,$
            vega=(keyword_set(abmag) eq 0),band_shift=0.0,/silent,_extra=extra)

          measure[sortindx[gthese[good]]].cflux = cflux
          measure[sortindx[gthese[good]]].uvflux = uvflux

          measure[sortindx[gthese[good]]].ubvrijhk_absmag = absmag
          measure[sortindx[gthese[good]]].ubvrijhk_absmag_ivar = absmag_ivar
          measure[sortindx[gthese[good]]].ubvrijhk_kcorrect = kcorr

; FUV,NUV, band_shift=0.0, AB
          kcorr = im_simple_kcorrect(measure[sortindx[gthese[good]]].zobj,$
            measure[sortindx[gthese[good]]].maggies,measure[sortindx[gthese[good]]].ivarmaggies,$
            filterlist,galex_filterlist(),restwave,restflux,absmag=absmag,$
            ivarabsmag=absmag_ivar,vega=0,band_shift=0.0,/silent,_extra=extra)

          measure[sortindx[gthese[good]]].galex_absmag = absmag
          measure[sortindx[gthese[good]]].galex_absmag_ivar = absmag_ivar
          measure[sortindx[gthese[good]]].galex_kcorrect = kcorr
          
; ugriz, band_shift=0.1, AB
          kcorr = im_simple_kcorrect(measure[sortindx[gthese[good]]].zobj,$
            measure[sortindx[gthese[good]]].maggies,measure[sortindx[gthese[good]]].ivarmaggies,$
            filterlist,ugriz_filterlist,restwave,restflux,absmag=absmag,$
            ivarabsmag=absmag_ivar,vega=0,band_shift=0.1,/silent,_extra=extra)

          measure[sortindx[gthese[good]]].ugriz_absmag = absmag
          measure[sortindx[gthese[good]]].ugriz_absmag_ivar = absmag_ivar
          measure[sortindx[gthese[good]]].ugriz_kcorrect = kcorr
       endif
       if (gchunk eq 0) then splog, format='("First '+$
         'Chunk = ",G0," minutes          ")', (systime(1)-t1)/60.0
    endfor
    splog, format='("All Chunks = ",G0," minutes          ")', $
      (systime(1)-t0)/60.0

; write out
    if (keyword_set(nowrite) eq 0)  then $
      im_mwrfits, measure, measurefile, /clobber
    
return
end
