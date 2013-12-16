;+
; NAME:
;   ISEDFIT_CHI2()
; PURPOSE:
;   Compute chi^2 and the best-fitting scale factor given a grid of
;   observed and synthetic photometry.  This routine is an ISEDFIT
;   support routine and in general should not be called on its own. 
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 10, NYU - developed
;   jm13aug09siena - updated to the latest data model 
;-

function isedfit_chi2, maggies, ivarmaggies, modelmaggies, $
  modelage=modelage, maxage=maxage, zindx=zindx, gchunk=gchunk, $
  ngalchunk=ngalchunk, ichunk=ichunk, nchunk=nchunk, $
  nminphot=nminphot, nzz=nzz, allages=allages, maxold=maxold, $
  silent=silent

    ndim = size(maggies,/n_dim)
    dims = size(maggies,/dim)
    nfilt = dims[0] ; number of filters
    if (ndim eq 1) then ngal = 1 else ngal = dims[1] ; number of galaxies

    nmodeldim = size(modelmaggies,/n_dim)
    modeldims = size(modelmaggies,/dim)
    if nmodeldim eq 3 then nmodel = modeldims[2] else nmodel = modeldims[1]
    
    gridchunk = {totalmass: -1.0, totalmass_err: -1.0, $
      chi2: 1E6, bestmaggies: fltarr(nfilt)}
    gridchunk = replicate(gridchunk,nmodel,ngal)

;   t0 = systime(1)
    for igal = 0L, ngal-1 do begin
       if (keyword_set(silent) eq 0) then begin
          if ((igal mod 50) eq 0) then print, format='("GalaxyChunk ",I0,"/",I0,", '+$
            'Modelgrid ",I0,"/",I0,", Galaxy ",I0,"/",I0,A10,$)', gchunk+1, ngalchunk, $
            ichunk+1, nchunk, igal, ngal, string(13b)
;         if ((igal mod 50) eq 0) then print, format='("GalaxyChunk ",I0,"/",I0,", '+$
;           'Modelgrid ",I0,"/",I0,", Galaxy ",I0,"/",I0,A10,$)', gchunk+1, ngalchunk, $
;           ichunk+1, nchunk, igal, ngal, string(13b)
       endif
; require detections or upper limits (i.e., non-zero ivar) in at least
; NMINPHOT bandpasses, and a non-zero flux and ivar in at least one
; bandpass, to set the overall normalization of the SED
       if (total(ivarmaggies[*,igal] gt 0.0) lt nminphot) or $
         (total((maggies[*,igal] gt 0.0) and (ivarmaggies[*,igal] gt 0.0)) eq 0.0) then continue
       nmaggies = 1D*abs(maggies[*,igal]) ; need absolute value to deal with negative fluxes correctly
       nivarmaggies = 1D*ivarmaggies[*,igal]
       dof = total(nivarmaggies gt 0)-1.0 ; degrees of freedom
       if (dof le 0) then message, 'This should not happen!'
;      t1 = systime(1)

       if keyword_set(maxold) then begin
          agediff = min(abs(modelage-maxage[igal]),these)
          nthese = 1
       endif else begin
; constrain the models to be younger than the universe
          if keyword_set(allages) then begin
             nthese = nmodel
             these = lindgen(nthese)
          endif else begin
             these = where((modelage le maxage[igal]),nthese)
          endelse
       endelse 
       if (nthese ne 0L) then begin ; at least one model
; interpolate the model photometry at the galaxy redshift
;         plot, modelgrid[imodel].modelmaggies[5,these,*],zindx[igal])*1.0D
          if nzz eq 1 then begin ; special case
             modelmaggies1 = interpolate(modelmaggies[*,these],findgen(nthese))
          endif else begin
             modelmaggies1 = interpolate(modelmaggies[*,*,these],$
               findgen(nfilt),zindx[igal],findgen(nthese),/grid)
          endelse
; perform acrobatic dimensional juggling to get the maximum likelihood
; scale-factor (total mass) and corresponding chi2 as a function of
; age; VSCALE is the maximum likelihood value of TOTALMASS and
; VSCALE_ERR is the 1-sigma error (see pg 84 of
; http://www.hep.phy.cam.ac.uk/~thomson/lectures/statistics/FittingHandout.pdf)
          vmodelmaggies = reform(1D*modelmaggies1,nfilt,nthese)
          vmaggies = rebin(reform(nmaggies,nfilt,1),nfilt,nthese)
          vivarmaggies = rebin(reform(nivarmaggies,nfilt,1),nfilt,nthese)
          vscale = total(reform((nivarmaggies*nmaggies),1,nfilt)#vmodelmaggies,1,/double)/$
            total(reform(nivarmaggies,1,nfilt)#vmodelmaggies^2,1,/double)
          vscale_err = 1.0/sqrt(total(reform(nivarmaggies,1,nfilt)#vmodelmaggies^2,1,/double))
          vchi2 = total(vivarmaggies*(vmaggies-rebin(reform(vscale,1,nthese),$
            nfilt,nthese)*vmodelmaggies)^2,1,/double)
          if total(vscale lt 0) ne 0.0 then message, 'Negative scale factor!'
; store the results
          bestmaggies = rebin(reform(vscale,1,nthese),nfilt,nthese)*vmodelmaggies
          gridchunk[these,igal].chi2 = vchi2/dof
          inf = where(finite(vscale) eq 0)
          if inf[0] ne -1 then message, 'Problem here'
          check = where((vscale le 0) or (finite(vscale) eq 0))
          if (check[0] ne -1) then message, 'Zero scale?!?'
          gridchunk[these,igal].totalmass = vscale            ; alog10(vscale)
          gridchunk[these,igal].totalmass_err = vscale_err    ; /vscale/alog(10)
          gridchunk[these,igal].bestmaggies = bestmaggies
       endif 
;      splog, format='("All models = ",G0," minutes")', (systime(1)-t1)/60.0       
    endfor                      ; galaxy loop
;   splog, format='("All galaxies = ",G0," minutes")', (systime(1)-t0)/60.0

return, gridchunk
end

