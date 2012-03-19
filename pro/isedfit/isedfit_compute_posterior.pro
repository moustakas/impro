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

function isedfit_compute_posterior, isedfit, modelgrid, fullgrid, $
  isedfit_post=isedfit_post

    ngal = n_elements(isedfit)
    nmodel = n_elements(modelgrid)
    nage = n_elements(modelgrid[0].age)
    nallmodel = nmodel*nage
;   nmodels = cmproduct((size(fullgrid,/dim))[0:1])

    ndraw = isedfit_ndraw() ; number of random draws
    allindx = lindgen(nallmodel)

; build the model grid parameter arrays we care about
    bigage    = reform(modelgrid.age,nallmodel)
    bigmass   = reform(modelgrid.mstar,nallmodel)
    bigtau    = reform(rebin(reform(modelgrid.tau,1,nmodel),nage,nmodel),nallmodel)
    bigZ      = reform(rebin(reform(modelgrid.Z,1,nmodel),nage,nmodel),nallmodel)
    bigmu     = reform(rebin(reform(modelgrid.mu,1,nmodel),nage,nmodel),nallmodel)
    bigav     = reform(rebin(reform(modelgrid.av,1,nmodel),nage,nmodel),nallmodel)
;   bigav     = reform(rebin(reform(modelgrid.mu*modelgrid.av,1,nmodel),nage,nmodel),nallmodel)

;   bigfburst   = reform(rebin(reform(modelgrid.fburst,1,nmodel),nage,nmodel),nallmodel)
;   bigtburst   = reform(rebin(reform(modelgrid.tburst,1,nmodel),nage,nmodel),nallmodel)
;   bigdtburst   = reform(rebin(reform(modelgrid.dtburst,1,nmodel),nage,nmodel),nallmodel)
;   bignburst   = reform(rebin(reform(modelgrid.nburst,1,nmodel),nage,nmodel),nallmodel)
;   bigtautrunc = reform(rebin(reform(modelgrid.tautrunc,1,nmodel),nage,nmodel),nallmodel)
    
; reconstruct the SFH, time-averaged SFR, and birthrate parameter for
; each model 
    bigsfr = bigage*0D
    bigsfr100 = bigage*0D ; average over the previous 100 Myr
    bigsfrage = bigage*0D ; SFR-weighted age
    bigb100 = bigage*0D   ; birthrate parameter
    bigmgal = bigage*0D   ; galaxy mass ignoring mass loss 
    for imod = 0L, nmodel-1 do begin
       tindx = lindgen(nage)+imod*nage
       sfr = isedfit_reconstruct_sfh(modelgrid[imod],outage=bigage[tindx],$
         sfr100=sfr100,sfrage=sfrage,b100=b100,mgalaxy=mgal)
       bigsfr[tindx] = sfr>1D-30 ; avoid SFR=0, which leads to an inf when we take the log, below 
       bigsfr100[tindx] = sfr100>1D-30
       bigb100[tindx] = b100
       bigmgal[tindx] = mgal
       bigsfrage[tindx] = sfrage
    endfor

; gotta loop...    
    for igal = 0L, ngal-1 do begin
       galgrid = reform(fullgrid[*,*,igal],nallmodel)

; get Monte Carlo draws from the posterior distributions of each
; parameter
       if (min(galgrid.chi2) lt 0.9D6) then begin
; additional priors
          prior = bigsfr*0+1
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
          isedfit_post[igal].scale = galgrid[allow[these]].scale
          isedfit_post[igal].scale_err = galgrid[allow[these]].scale_err
          isedfit_post[igal].chi2 = galgrid[allow[these]].chi2

; multiply the stellar masses and SFRs of the models by the scale
; factor; perturb the scale factor by the Gaussian error to avoid the
; posterior distribution being infinitely sharp
          logscale_err = galgrid[allow[these]].scale_err/galgrid[allow[these]].scale/alog(10)
          logscale = alog10(galgrid[allow[these]].scale)+randomn(seed,ndraw)*logscale_err

          isedfit[igal] = isedfit_packit(isedfit[igal],alog10(bigmass[allow[these]])+logscale,type='mass')
          isedfit[igal] = isedfit_packit(isedfit[igal],alog10(bigsfr[allow[these]])+logscale,type='sfr')
          isedfit[igal] = isedfit_packit(isedfit[igal],alog10(bigsfr100[allow[these]])+logscale,type='sfr100')

          isedfit[igal] = isedfit_packit(isedfit[igal],bigsfrage[allow[these]],type='sfrage')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigage[allow[these]],type='age')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigtau[allow[these]],type='tau')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigZ[allow[these]],type='Z')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigav[allow[these]],type='av')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigmu[allow[these]],type='mu')
          isedfit[igal] = isedfit_packit(isedfit[igal],bigb100[allow[these]],type='b100')

          neg = where(isedfit_post[igal].scale le 0)
          if (neg[0] ne -1) then message, 'Negative scale factor!'
       
; now compute the maximum likelihood (best-fit) estimates of the
; various parameters
          ml = max(post,mindx)
          isedfit[igal].chi2 = galgrid[mindx].chi2
          isedfit[igal].scale = galgrid[mindx].scale
          isedfit[igal].scale_err = galgrid[mindx].scale_err
          
          allmindx = array_indices([nage,nmodel],mindx,/dim) ; parse the index
          isedfit[igal].ageindx = allmindx[0]
          isedfit[igal] = im_struct_assign(modelgrid[allmindx[1]],isedfit[igal],/nozero)

;         isedfit[igal].modelindx = modelgrid[allmindx[1]].modelindx
;         isedfit[igal].chunkindx = modelgrid[allmindx[1]].chunkindx
;         isedfit[igal].delayed = modelgrid[allmindx[1]].delayed
;         isedfit[igal].bursttype = modelgrid[allmindx[1]].bursttype
;         
;         isedfit[igal].nburst = modelgrid[allmindx[1]].nburst
;         isedfit[igal].tautrunc = modelgrid[allmindx[1]].tautrunc
;         isedfit[igal].tburst = modelgrid[allmindx[1]].tburst
;         isedfit[igal].dtburst = modelgrid[allmindx[1]].dtburst
;         isedfit[igal].fburst = modelgrid[allmindx[1]].fburst
          
;         isedfit[igal].tau = bigtau[mindx]
;         isedfit[igal].Z = bigZ[mindx]
;         isedfit[igal].av = bigav[mindx]
;         isedfit[igal].mu = bigmu[mindx]

          isedfit[igal].age = bigage[mindx]
          isedfit[igal].sfrage = bigsfrage[mindx]
          isedfit[igal].mass = alog10(bigmass[mindx]*isedfit[igal].scale)
          isedfit[igal].sfr = alog10(bigsfr[mindx]*isedfit[igal].scale)
          isedfit[igal].sfr100 = alog10(bigsfr100[mindx]*isedfit[igal].scale)
          isedfit[igal].b100 = bigb100[mindx]
          isedfit[igal].bestmaggies = galgrid[mindx].bestmaggies

;         psfile = 'chi2_vs_deltamag.ps'
;         im_plotconfig, 0, pos, psfile=psfile, height=5.5
;         dmag = 2.5*alog10(galgrid.bestmaggies[17])-2.5*alog10(isedfit.maggies[17])
;         djs_plot, dmag, galgrid.chi2, psym=6, xtitle='[ch2]_{model}-[ch2]_{data}', $
;           ytitle='\chi^{2}_{\nu}', symsize=0.5, position=pos, xsty=1, ysty=1, thick=1, $
;           yrange=[0,6], color=im_color('grey')
;         djs_oplot, dmag[allow[these]], (galgrid.chi2)[allow[these]], psym=6, $
;           symsize=0.5, thick=1, color=im_color('firebrick')
;         im_plotconfig, psfile=psfile, /psclose, /pdf
       endif
; some plots
;       xx = bigage & yy = bigdtburst
;       djs_plot, xx, galgrid.chi2, ps=6, /ylog, xlog=0, yr=isedfit[igal].chi2*[0.9,3], xsty=3, ysty=3, sym=0.5 & djs_oplot, xx[allow[these]], galgrid[allow[these]].chi2, ps=6, sym=0.5, color='green'
;
;      ww = where(galgrid.chi2 lt 150) & w2 = where(galgrid.chi2 lt 100)
;      djs_plot, xx, yy, ps=6, /ylog, sym=0.2, ysty=3 & djs_oplot, xx[ww], yy[ww], psym=7, color='red', sym=2 & djs_oplot, xx[w2], yy[w2], psym=7, color='green', sym=2.5
    endfor 

return, isedfit
end
