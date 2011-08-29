;+
; NAME:
;   ISEDFIT_RESTORE()
;
; PURPOSE:
;   Function to restore the best-fitting iSEDfit model. 
;
; INPUTS:
;   paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   iopath - I/O path
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

function isedfit_restore, paramfile, isedfit, params=params, iopath=iopath, $
  index=index, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix, $
  flambda=flambda, fnu=fnu, nomodels=nomodels, silent=silent

    if (n_elements(paramfile) eq 0L) and (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_restore'
       return, -1
    endif

    if (n_elements(iopath) eq 0) then iopath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)
    fp = isedfit_filepaths(params,outprefix=outprefix,iopath=iopath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

; restore the ISEDFIT output
    if (file_test(fp.iopath+fp.isedfit_outfile+'.gz',/regular) eq 0L) then $
      message, 'ISEDFIT output not found!'
    splog, 'Reading '+fp.iopath+fp.isedfit_outfile+'.gz'
    isedfit1 = mrdfits(fp.iopath+fp.isedfit_outfile+'.gz',1,/silent)

; only restore a subset of the objects    
    if (n_elements(index) eq 0L) then isedfit = isedfit1 else $
      isedfit = isedfit1[index]
    ngal = n_elements(isedfit)

    if keyword_set(nomodels) then begin
       return, 1
    endif else begin
       light = 2.99792458D18             ; speed of light [A/s]
       dist = 10.0*3.085678D18           ; fiducial distance [10 pc in cm]
       dlum = dluminosity(isedfit.zobj,/cm) ; luminosity distance [cm]

; read the models; group by chunks; remove objects that were not
; fitted
       allchunks = isedfit.chunkindx
;      crap = where(allchunks eq -1,ncrap)
;      if (ncrap ne 0) then allchunks[crap] = 0 ; HACK!!
       chunks = allchunks[uniq(allchunks,sort(allchunks))]
       nchunk = n_elements(chunks)
       
; initialize the model structure (need to get the number of pixels)
       junk = mrdfits(fp.sfhgrid_chunkfiles[0],1,row=0,/silent)
       npix = n_elements(junk.wave)

       model = {wave: fltarr(npix), flux: fltarr(npix)}
       model = replicate(temporary(model),ngal)
       
       for ichunk = 0L, nchunk-1L do begin
          these = where(chunks[ichunk] eq allchunks,nthese)
          if (nthese ne 0L) and (chunks[ichunk] ge 0) then begin
             chunkfile = strtrim(fp.sfhgrid_chunkfiles[chunks[ichunk]],2)
             if (not keyword_set(silent)) then $
               splog, 'Reading '+chunkfile
             grid = mrdfits(chunkfile,1,/silent)
             for ii = 0L, nthese-1L do begin
                if (isedfit[these[ii]].chi2 lt 1E6) then begin
                   i1 = isedfit[these[ii]].modelindx; mod 100
                   i2 = isedfit[these[ii]].ageindx
                   model[these[ii]].wave = grid[i1].wave
                   model[these[ii]].flux = grid[i1].flux[*,i2]
                endif
             endfor
          endif  
       endfor

; read the IGM absorption table
       if params.igm then begin
          splog, 'Reading IGM attenuation lookup table'
          igmgrid = mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)
       endif

; loop through every object and construct the best-fit model
       for igal = 0L, ngal-1L do begin
          if (isedfit[igal].chi2 lt 1E6) then begin
             zobj = isedfit[igal].zobj
             zwave = model[igal].wave*(1+zobj)
             zflux_flam = isedfit[igal].scale*model[igal].flux*$ ; [erg/s/cm2/A]
               (dist/dlum[igal])^2.0/(1.0+zobj)
             if params.igm then begin
                windx = findex(igmgrid.wave,zwave)
                zindx = findex(igmgrid.zgrid,zobj)
                igm = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
                zflux_flam = zflux_flam*igm
             endif
             zflux_fnu = zflux_flam*zwave^2/light ; [erg/s/cm2/Hz]
             zflux_ab = -2.5*alog10(zflux_fnu>1D-50)-48.6
             model[igal].wave = zwave
             model[igal].flux = zflux_ab ; default
             if keyword_set(fnu) then model[igal].flux = zflux_fnu
             if keyword_set(flambda) then model[igal].flux = zflux_flam
          endif
       endfor
    endelse

return, model
end
