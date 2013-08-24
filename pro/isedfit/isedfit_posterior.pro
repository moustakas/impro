;+
; NAME:
;   ISEDFIT_POSTERIOR()
; PURPOSE:
;   Compute the posterior probability distribution.  This routine is
;   an ISEDFIT support routine and in general should not be called on
;   its own.  
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 10, NYU - developed
;-

function isedfit_packit, isedfit, array, type=type, wquant=wquant
; support routine for isedfit_posterior()
    
    tag50 = tag_indx(isedfit,type+'_50')
    tagerr = tag_indx(isedfit,type+'_err')
    tagavg = tag_indx(isedfit,type+'_avg')
;   tagefferr = tag_indx(isedfit,type+'_eff_err')

    isedfit.(tagavg) = djs_mean(array)
;   isedfit.(tagerr) = djsig(array)

    quant = [1.0-gauss_pdf(2.0),0.5,gauss_pdf(2.0)]
    wquant = weighted_quantile(array,quant=quant)
    isedfit.(tag50) = wquant[1]
    isedfit.(tagerr) = (wquant[2]-wquant[0])/4.0
;   isedfit.(tagefferr) = (wquant[2]-wquant[0])/4.0

return, isedfit
end

function isedfit_posterior, isedfit, modelgrid=modelgrid, $
  galaxygrid=galaxygrid, params=params, isedfit_post=isedfit_post

    ngal = n_elements(isedfit)
    nmodel = n_elements(modelgrid)
    ndraw = params.ndraw ; number of random draws

; gotta loop...    
    for igal = 0L, ngal-1 do begin
       if ((igal mod 50) eq 0) then print, format='("Computing the posterior '+$
         'for galaxy ",I0,"/",I0, A10,$)', igal, ngal, string(13b)

; get Monte Carlo draws from the posterior distributions of each
; parameter
       if (min(galaxygrid[*,igal].chi2) lt 0.9D6) then begin
; could add additional priors here
          prior = modelgrid.sfr*0+1
          post = prior*exp(-0.5D*(galaxygrid[*,igal].chi2-min(galaxygrid[*,igal].chi2)))
          if (total(post,/double) eq 0.0) then begin
             splog, 'Chi^2 value too large for object '+strtrim(igal,2)
             continue
          endif
          post = post/total(post,/double)
          allow = where(galaxygrid[*,igal].chi2 lt 1E6,nallow) ; ignore rejected models

; need to pad with zero probability because we're using LONG()
          these = long(interpolate(lindgen(nallow+1),findex([0D,$
            total(post[allow],/cumu,/double)],randomu(seed,ndraw))))

          isedfit_post[igal].draws = allow[these]
          isedfit_post[igal].totalmass = galaxygrid[allow[these],igal].totalmass
          isedfit_post[igal].totalmass_err = galaxygrid[allow[these],igal].totalmass_err
          isedfit_post[igal].chi2 = galaxygrid[allow[these],igal].chi2

; multiply the stellar masses and SFRs of the models by the scale
; factor (total mass); perturb the scale factor by the Gaussian error
; to avoid the posterior distribution being infinitely sharp
          logscale_err = galaxygrid[allow[these],igal].totalmass_err/galaxygrid[allow[these],igal].totalmass/alog(10)
          logscale = alog10(galaxygrid[allow[these],igal].totalmass)+randomn(seed,ndraw)*logscale_err

          isedfit[igal] = isedfit_packit(isedfit[igal],alog10(modelgrid[allow[these]].mstar)+logscale,type='mstar')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].sfr+logscale,type='sfr')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].sfr100+logscale,type='sfr100')
          
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].b100,type='b100')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].age,type='age')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].sfrage,type='sfrage')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].tau,type='tau')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].Zmetal,type='Zmetal')

          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].av,type='av')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].mu,type='mu')

          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].oiiihb,type='oiiihb')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].ewoii,type='ewoii')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].ewoiiihb,type='ewoiiihb')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].ewniiha,type='ewniiha')
          isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow[these]].oiiihb,type='oiiihb')
          
          neg = where(isedfit_post[igal].totalmass le 0)
          if (neg[0] ne -1) then message, 'Negative totalmass factor!'

; now compute the maximum likelihood (best-fit) estimates of the
; various parameters
          ml = max(post,mindx)
          isedfit[igal].chi2 = galaxygrid[mindx,igal].chi2
          isedfit[igal].totalmass = galaxygrid[mindx,igal].totalmass
          isedfit[igal].totalmass_err = galaxygrid[mindx,igal].totalmass_err
          isedfit[igal].bestmaggies = galaxygrid[mindx,igal].bestmaggies
          isedfit[igal] = im_struct_assign(modelgrid[mindx],isedfit[igal],/nozero)

          isedfit[igal].age = modelgrid[mindx].age
          isedfit[igal].sfrage = modelgrid[mindx].sfrage
          isedfit[igal].mstar = alog10(modelgrid[mindx].mstar*isedfit[igal].totalmass)
          isedfit[igal].sfr = modelgrid[mindx].sfr+alog10(isedfit[igal].totalmass)
          isedfit[igal].sfr100 = modelgrid[mindx].sfr100+alog10(isedfit[igal].totalmass)
          isedfit[igal].b100 = modelgrid[mindx].b100
       endif
    endfor 

return, isedfit
end
