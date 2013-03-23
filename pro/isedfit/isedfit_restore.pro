;+
; NAME:
;   ISEDFIT_RESTORE()
;
; PURPOSE:
;   Restore the best-fitting iSEDfit model. 
;
; INPUTS:
;   paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - iSEDfit parameter data structure (over-rides PARAMFILE) 
;   isedpath - I/O path
;   index - zero-indexed list of objects to restore (default is to
;     restore all)
;   isedfit_sfhgrid_dir - see BUILD_ISEDFIT_SFHGRID
;   outprefix - optional output prefix string (see ISEDFIT) 
;   in_isedfit - like ISEDFIT, but as an input; needed to enable some
;     of the fancy plotting I've done for papers
;
; KEYWORD PARAMETERS:
;   flambda - convert to F(lambda) (erg/s/cm^2/A) (default is AB mag) 
;   fnu - convert to F(nu) (erg/s/cm^2/Hz) (default is AB mag) 
;   nomodels - restore the ISEDFIT structure, but not the model
;     spectra  
;   noigm - do not convolve with the IGM
;   silent - suppress messages to STDOUT
;
; OUTPUTS:
;   model - output data structure array [NGAL]
;     .WAVE - observed-frame wavelength array (Angstrom)
;     .FLUX - observed-frame flux array (AB mag)
;
; OPTIONAL OUTPUTS:
;   isedfit - iSEDfit result structure (see ISEDFIT) [NGAL] 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 27, NYU - largely excised from
;     ISEDFIT_QAPLOT and ISEDFIT_MEASURE
;   jm11oct01ucsd - documentation cleaned up and updated
;
; Copyright (C) 2007, 2011, John Moustakas
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

function isedfit_restore, paramfile, isedfit, params=params, isedfit_dir=isedfit_dir, $
  supergrid_paramfile=supergrid_paramfile, thissupergrid=thissupergrid, $
  montegrids_dir=montegrids_dir, index=index, outprefix=outprefix, flambda=flambda, $
  fnu=fnu, nomodels=nomodels, noigm=noigm, silent=silent, in_isedfit=in_isedfit

    if (n_elements(paramfile) eq 0L) and (n_elements(params) eq 0) then begin
       doc_library, 'isedfit_restore'
       return, -1
    endif

    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile)
    fp = isedfit_filepaths(params,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=thissupergrid,isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir,$
      outprefix=outprefix)

; restore the ISEDFIT output; optionally take ISEDFIT as an input
; structure
    if (n_elements(in_isedfit) eq 0L) then begin
       if (file_test(fp.isedfit_dir+fp.isedfit_outfile+'.gz',/regular) eq 0L) then $
         message, 'ISEDFIT output not found!'
       splog, 'Reading '+fp.isedfit_dir+fp.isedfit_outfile+'.gz'
       isedfit1 = mrdfits(fp.isedfit_dir+fp.isedfit_outfile+'.gz',1,/silent)

; only restore a subset of the objects    
       if (n_elements(index) eq 0L) then isedfit = isedfit1 else $
         isedfit = isedfit1[index]
    endif else isedfit = in_isedfit
    ngal = n_elements(isedfit)

    if keyword_set(nomodels) then begin
       return, 1
    endif else begin
       light = 2.99792458D18             ; speed of light [A/s]
       pc10 = 3.085678D19                ; fiducial distance [10 pc in cm]
       dlum = pc10*10D^(lf_distmod(isedfit.zobj,omega0=params.omega0,$ ; [cm]
         omegal0=params.omegal)/5D)/params.h100 
;      dlum = dluminosity(isedfit.zobj,/cm) ; luminosity distance [cm]
       
; read the models; group by chunks; remove objects that were not
; fitted
       allchunks = isedfit.chunkindx
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
       if params.igm and (keyword_set(noigm) eq 0) then begin
          splog, 'Reading IGM attenuation lookup table'
          igmgrid = mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)
       endif

; loop through every object and construct the best-fit model
       for igal = 0L, ngal-1L do begin
          if (isedfit[igal].chi2 lt 1E6) then begin
             zobj = isedfit[igal].zobj
             zwave = model[igal].wave*(1+zobj)
             zflux_flam = isedfit[igal].scale*model[igal].flux*$ ; [erg/s/cm2/A]
               (pc10/dlum[igal])^2.0/(1.0+zobj)
             if params.igm and (keyword_set(noigm) eq 0) then begin
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
