;+
; NAME:
;   ISEDFIT_COMPUTE_CHI2()
;
; PURPOSE:
;   Compute chi^2 and the best-fitting scale factor (stellar mass)
;   given a grid of observed and synthetic photometry.  This routine
;   is an ISEDFIT support routine, and in general should not be called
;   on its own.
;
; INPUTS: 
;   maggies - input photometry [NFILT,NGAL]
;   ivarmaggies - corresponding inverse variances [NFILT,NGAL] 
;   chunkmodels - see ISEDFIT
;   maxage - age of the universe at the redshift of the galaxy [NGAL] 
;   zindx - see ISEDFIT
;
; OPTIONAL INPUTS: 
;   The following variables are passed from ISEDFIT so that an
;   informational message can be printed to STDOUT:
;     gchunk
;     ngalchunk
;     ichunk
;     nchunk
;   nminphot - an object must have good photometry in at least
;     NMINPHOT bandpasses
;
; KEYWORD PARAMETERS: 
;   allages - by default only consider models that are younger than
;     the age of the universe at the redshift of the galaxy, unless
;     /ALLAGES 
;
; OUTPUTS: 
;   gridchunk - output data structure containing the chi^2, stellar
;     mass, and best-fitting synthesized photometry for each
;     combination of galaxy and models [NMODEL,NAGE,NGAL]
;
; OPTIONAL OUTPUTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb - developed
;
; Copyright (C) 2009, John Moustakas
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

function isedfit_compute_chi2, maggies, ivarmaggies, chunkmodels, maxage, $
  zindx, gchunk=gchunk, ngalchunk=ngalchunk, ichunk=ichunk, nchunk=nchunk, $
  nminphot=nminphot, allages=allages, maxold=maxold, silent=silent

    ndim = size(maggies,/n_dim)
    dims = size(maggies,/dim)
    nfilt = dims[0] ; number of filters
    if (ndim eq 1L) then ngal = 1L else ngal = dims[1] ; number of galaxies
    
    nmodel = n_elements(chunkmodels)
    nage = n_elements(chunkmodels[0].age)

    gridchunk = {mass: -1.0, mass_err: -1.0, chi2: 1E6, $ ; zptoffset: fltarr(nfilt), $
      bestmaggies: fltarr(nfilt)}
    gridchunk = replicate(gridchunk,nage,nmodel,ngal)

;   t0 = systime(1)
    for igal = 0L, ngal-1L do begin
       if (keyword_set(silent) eq 0) then begin
          if ((igal mod 50) eq 0) then print, format='("GalaxyChunk ",I0,"/",I0,", '+$
            'ModelChunk ",I0,"/",I0,", Galaxy ",I0,"/",I0,A10,$)', gchunk+1L, ngalchunk, $
            ichunk+1L, nchunk, igal+1L, ngal, string(13b)
       endif
; require detections or upper limits (i.e., non-zero ivar) in at least
; NMINPHOT bandpasses, and a non-zero flux and ivar in at least one
; bandpass, to set the overall normalization of the SED
       if (total(ivarmaggies[*,igal] gt 0.0) lt nminphot) or $
         (total((maggies[*,igal] gt 0.0) and (ivarmaggies[*,igal] gt 0.0)) eq 0.0) then continue
       nmaggies = abs(maggies[*,igal]*1.0D) ; need absolute value to deal with negative fluxes correctly
       nivarmaggies = ivarmaggies[*,igal]*1.0D
;      t1 = systime(1)
       for imodel = 0L, nmodel-1L do begin
          if keyword_set(maxold) then begin
             agediff = min(abs(chunkmodels[imodel].age-maxage[igal]),these)
             nthese = 1
          endif else begin
; constrain the models to be younger than the universe
             if keyword_set(allages) then begin
                nthese = nage
                these = lindgen(nthese)
             endif else begin
                these = where((chunkmodels[imodel].age le maxage[igal]),nthese)
             endelse
          endelse
          if (nthese eq 0) then message, 'Your galaxy is too young!'
; interpolate the model photometry at the galaxy redshift
;         plot, chunkmodels[imodel].modelmaggies[5,these,*],zindx[igal])*1.0D
          modelmaggies = interpolate(chunkmodels[imodel].modelmaggies[*,these,*],zindx[igal])*1.0D
; perform acrobatic dimensional juggling to get the maximum likelihood
; mass and corresponding chi2 as a function of age; VMASS is the
; maximum likelihood mass and VMASSERR is the 1-sigma error (see pg 84
; of http://www.hep.phy.cam.ac.uk/~thomson/lectures/statistics/FittingHandout.pdf)
          vmodelmaggies = reform(modelmaggies,nfilt,nthese)
          vmaggies = rebin(reform(nmaggies,nfilt,1),nfilt,nthese)
          vivarmaggies = rebin(reform(nivarmaggies,nfilt,1),nfilt,nthese)
          vmass = total(reform((nivarmaggies*nmaggies),1,nfilt)#vmodelmaggies,1,/double)/$
            total(reform(nivarmaggies,1,nfilt)#vmodelmaggies^2,1,/double)
          vmass_err = 1.0/sqrt(total(reform(nivarmaggies,1,nfilt)#vmodelmaggies^2,1,/double))
          vchi2 = total(vivarmaggies*(vmaggies-rebin(reform(vmass,1,nthese),$
            nfilt,nthese)*vmodelmaggies)^2,1,/double)
;         print, total(nivarmaggies*(nmaggies-(vmass[30]+vmass_err[30])*vmodelmaggies[*,30])^2)-vchi2[30]

;; for some galaxy samples maggies can be formally negative, and
;; therefore the normalization (i.e., the mass) can be negative too; if
;; this happens then issue a warning message and take the absolute
;; value
;          neg = where(vmass lt 0.0,nneg)
;          if (nneg ne 0L) then begin
;             if (imodel eq 0L) then splog, 'Formally negative '+$
;               'stellar mass for galaxy index '+strtrim(igal,2)
;             vmass = abs(vmass)
;          endif
; store the results
          bestmaggies = rebin(reform(vmass,1,nthese),nfilt,nthese)*vmodelmaggies
          gridchunk[these,imodel,igal].chi2 = vchi2
          inf = where(finite(vmass) eq 0)
          if inf[0] ne -1 then message, 'Problem here'
          check = where((vmass le 0) or (finite(vmass) eq 0))
          if (check[0] ne -1) then message, 'Zero mass?!?'
          gridchunk[these,imodel,igal].mass = alog10(vmass)
          gridchunk[these,imodel,igal].mass_err = vmass_err/vmass/alog(10)
          gridchunk[these,imodel,igal].bestmaggies = bestmaggies
       endfor                   ; model loop
;      splog, format='("All models = ",G0," minutes")', (systime(1)-t1)/60.0       
    endfor                      ; galaxy loop
;   splog, format='("All galaxies = ",G0," minutes")', (systime(1)-t0)/60.0

return, gridchunk
end

