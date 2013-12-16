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

function isedfit_packit, isedfit, array, weight=weight, $
  type=type, wquant=wquant
; support routine for isedfit_posterior()
    
    tag50 = tag_indx(isedfit,type+'_50')
    tagerr = tag_indx(isedfit,type+'_err')
    tagavg = tag_indx(isedfit,type+'_avg')
;   tagefferr = tag_indx(isedfit,type+'_eff_err')
    
    isedfit.(tagavg) = im_weighted_mean(array,weight=weight)
;   isedfit.(tagavg) = djs_mean(array)
;   isedfit.(tagerr) = djsig(array)

    quant = [1.0-gauss_pdf(2.0),0.5,gauss_pdf(2.0)]
    wquant = weighted_quantile(array,weight,quant=quant)
    isedfit.(tag50) = wquant[1]
    isedfit.(tagerr) = (wquant[2]-wquant[0])/4.0
;   isedfit.(tagefferr) = (wquant[2]-wquant[0])/4.0

return, isedfit
end

function isedfit_posterior, isedfit, params=params, modelgrid=modelgrid, $
  galaxygrid=galaxygrid, isedfit_post=isedfit_post, nmin=nmin, $
  bestfitonly=bestfitonly, silent=silent

    ngal = n_elements(isedfit)
    nmodel = n_elements(modelgrid)
    ndraw = params.ndraw ; number of random draws

    pars = ['b100','b1000','age','sfrage','tau','Zmetal',$
      'av','mu','oiiihb','ewoii','ewoiiihb','ewniiha']
    npars = n_elements(pars)
    
; gotta loop...    
    for igal = 0L, ngal-1 do begin
       if keyword_set(silent) eq 0 then if ((igal mod 50) eq 0) then $
         print, format='("Computing the posterior for galaxy ",I0,"/",I0, A10,$)', $
         igal, ngal, string(13b)

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
          postnorm = total(post,/double)
          post = post/postnorm
          allow = where(galaxygrid[*,igal].chi2 lt 1E6,nallow) ; ignore rejected models

; first get the maximum likelihood (best-fit) estimates
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
          isedfit[igal].b1000 = modelgrid[mindx].b1000

; next get the moments of the posterior; for the stellar masses and
; SFRs we have to multiply the models by the scale factor (total
; mass); perturb the scale factor by the Gaussian error to avoid the
; posterior distribution being infinitely sharp
          if keyword_set(bestfitonly) eq 0 then begin
             logscale_err = galaxygrid[allow,igal].totalmass_err/$
               galaxygrid[allow,igal].totalmass/alog(10)
             logscale = alog10(galaxygrid[allow,igal].totalmass)+$
               randomn(seed,nallow)*logscale_err

             isedfit[igal] = isedfit_packit(isedfit[igal],alog10(modelgrid[allow].mstar)+$
               logscale,weight=post[allow],type='mstar')
             isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow].sfr+$
               logscale,weight=post[allow],type='sfr')
             isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow].sfr100+$
               logscale,weight=post[allow],type='sfr100')

             for ip = 0, npars-1 do begin
                pindx = tag_indx(modelgrid[0],pars[ip])
                isedfit[igal] = isedfit_packit(isedfit[igal],modelgrid[allow].(pindx),$
                  weight=post[allow],type=pars[ip])
             endfor
             
; randomly sample the posterior for our output structure; need
; to pad with zero probability because we're using LONG()
             these = long(interpolate(lindgen(nallow+1),findex([0D,$
               total(post[allow],/cumu,/double)],randomu(seed,ndraw))))

             isedfit_post[igal].draws = allow[these]
             isedfit_post[igal].totalmass = galaxygrid[allow[these],igal].totalmass
             isedfit_post[igal].totalmass_err = galaxygrid[allow[these],igal].totalmass_err
             isedfit_post[igal].chi2 = galaxygrid[allow[these],igal].chi2
          endif 
       endif 
    endfor 

return, isedfit
end
