;+
; NAME:
;   ISEDFIT_MINIMIZE_CHI2()
;
; PURPOSE:
;   Internal support routine for ISEDFIT: minimize the chi^2 vector. 
;
; INPUTS: 
;   isedfit - 
;   modelgrid - 
;   fullgrid - 
;   nallmodel - 
;
; OPTIONAL INPUTS: 
;   quant - 
;
; KEYWORD PARAMETERS: 
;   debug - 
;
; OUTPUTS: 
;   isedfit - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2010, John Moustakas
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

function isedfit_packit, isedfit, array, type=type

    tagavg = tag_indx(isedfit,type+'_avg')
    tagerr = tag_indx(isedfit,type+'_err')
    tag50 = tag_indx(isedfit,type+'_50')
    tagefferr = tag_indx(isedfit,type+'_eff_err')

    isedfit.(tagavg) = djs_mean(array)
    isedfit.(tagerr) = djsig(array)

    quant = [1.0-gauss_pdf(2.0),0.5,gauss_pdf(2.0)]
    wquant = weighted_quantile(array,quant=quant)
    isedfit.(tag50) = wquant[1]
    isedfit.(tagefferr) = (wquant[2]-wquant[0])/4.0

return, isedfit
end

function isedfit_compute_posterior, isedfit, modelgrid, fullgrid, $
  isedfit_post=isedfit_post, debug=debug

    ngal = n_elements(isedfit)
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].age)
    nallmodel = nmodel*nage
;   nmodels = cmproduct((size(fullgrid,/dim))[0:1])

    ndraw = isedfit_ndraw() ; number of random draws
    allindx = lindgen(nallmodel)

; build the model grid parameter arrays we care about
    bigage    = reform(modelgrid.age,nallmodel)
    bigtau    = reform(rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel),nallmodel)
    bigZ      = reform(rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel),nallmodel)
    bigebv    = reform(rebin(reform(modelgrid.ebv,1,nmodel),nage,nmodel),nallmodel)
    bigmu     = reform(rebin(reform(modelgrid.mu,1,nmodel),nage,nmodel),nallmodel)
    bignburst = reform(rebin(reform(modelgrid.nburst,1,nmodel),nage,nmodel),nallmodel)
    
; reconstruct the SFH, time-averaged SFR, and birthrate parameter for
; each model 
    bigsfr = bigage*0D
    bigsfr100 = bigage*0D ; average over the previous 100 Myr
    bigb100 = bigage*0D   ; birthrate parameter
;   bigmgal = bigage*0D   ; galaxy mass ignoring mass loss 
    for imod = 0L, nmodel-1 do begin
       these = lindgen(nage)+imod*nage
       sfr = isedfit_reconstruct_sfh(modelgrid[imod],age=bigage[these],$
         sfr100=sfr100,b100=b100);,mgalaxy=mgal)
       bigsfr[these] = sfr
       bigsfr100[these] = sfr100
       bigb100[these] = b100
;      bigmgal[these] = mgal
    endfor

; additional priors: do not allow dusty, non-star-forming solutions
; (note that disallowing models that are too old is built into
; ISEDFIT_COMPUTE_CHI2) 
    prior = bigb100*0+1
;   prior = ((bigb100 lt 1D-2) and (bigebv gt 0.0)) eq 0

; gotta loop...    
    for igal = 0L, ngal-1 do begin
       galgrid = reform(fullgrid[*,*,igal],nallmodel)

; get Monte Carlo draws from the posterior distributions of each
; parameter
       post = prior*exp(-0.5D*(galgrid.chi2-min(galgrid.chi2)))
       if (total(post,/double) eq 0.0) then begin
          splog, 'Chi^2 value too large for object '+strtrim(igal,2)
          continue
       endif
       post = post/total(post,/double)
       allow = where(post gt 0.0,nallow)

; need to pad with zero probability because we're using LONG()
       these = long(interpolate(lindgen(nallow+1),findex([0D,$
         total(post[allow],/cumu,/double)],randomu(seed,ndraw))))
       isedfit_post[igal].draws = allow[these]
       isedfit_post[igal].mass = galgrid[allow[these]].mass+randomn(seed,ndraw)*$
         galgrid[allow[these]].mass_err

;; no need to store these distributions as we can easily reconstruct
;; them from the parent MONTEGRIDS       
;       isedfit_post[igal].Z      = bigZ[allow[these]]
;       isedfit_post[igal].tau    = bigtau[allow[these]]
;       isedfit_post[igal].age    = bigage[allow[these]]
;       isedfit_post[igal].ebv    = bigebv[allow[these]]
;       isedfit_post[igal].mu     = bigmu[allow[these]]
;       isedfit_post[igal].b100   = bigb100[allow[these]]
;       isedfit_post[igal].sfr    = 10.0^isedfit_post[igal].mass*bigsfr[allow[these]]
;       isedfit_post[igal].sfr100 = 10.0^isedfit_post[igal].mass*bigsfr100[allow[these]]

       isedfit[igal] = isedfit_packit(isedfit[igal],isedfit_post[igal].mass,type='mass')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigage[allow[these]],type='age')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigtau[allow[these]],type='tau')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigZ[allow[these]],type='Z')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigebv[allow[these]],type='ebv')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigmu[allow[these]],type='mu')
       isedfit[igal] = isedfit_packit(isedfit[igal],bigb100[allow[these]],type='b100')
       isedfit[igal] = isedfit_packit(isedfit[igal],10.0^isedfit_post[igal].mass*$
         bigsfr[allow[these]],type='sfr')
       isedfit[igal] = isedfit_packit(isedfit[igal],10.0^isedfit_post[igal].mass*$
         bigsfr100[allow[these]],type='sfr100')

       neg = where(isedfit_post[igal].mass le 0)
       if (neg[0] ne -1) then message, 'This is bad'
       
; now compute the maximum likelihood (best-fit) estimates of the
; various parameters
       ml = max(post,mindx)
       isedfit[igal].chi2 = galgrid[mindx].chi2

       allmindx = array_indices([nage,nmodel],mindx,/dim) ; parse the index
       isedfit[igal].ageindx = allmindx[0]
       isedfit[igal].imf = modelgrid[allmindx[1]].imf
       isedfit[igal].modelindx = modelgrid[allmindx[1]].modelindx
       isedfit[igal].chunkindx = modelgrid[allmindx[1]].chunkindx

       isedfit[igal].nburst = modelgrid[allmindx[1]].nburst
       isedfit[igal].tauburst = modelgrid[allmindx[1]].tauburst
       isedfit[igal].tburst = modelgrid[allmindx[1]].tburst
       isedfit[igal].dtburst = modelgrid[allmindx[1]].dtburst
       isedfit[igal].fburst = modelgrid[allmindx[1]].fburst

       isedfit[igal].tau = bigtau[mindx]
       isedfit[igal].Z = bigZ[mindx]
       isedfit[igal].ebv = bigebv[mindx]
       isedfit[igal].mu = bigmu[mindx]
       isedfit[igal].age = bigage[mindx]
       isedfit[igal].b100 = bigb100[mindx]

       isedfit[igal].mass = galgrid[mindx].mass
       isedfit[igal].bestmaggies = galgrid[mindx].bestmaggies

       isedfit[igal].sfr = 10.0^isedfit[igal].mass*bigsfr[mindx]
       isedfit[igal].sfr100 = 10.0^isedfit[igal].mass*bigsfr100[mindx]
    endfor 

return, isedfit
end
