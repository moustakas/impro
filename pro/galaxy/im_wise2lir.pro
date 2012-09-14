;+
; NAME:
;   IM_WISE2LIR()
;
; PURPOSE:
;   Given a redshift and 12- and 22-micron *observed* WISE photometry,
;   compute the total infrared luminosity, L(IR)=L(8-1100), based on
;   various infrared SED models.
;
; INPUTS: 
;   redshift - redshift for each object [NGAL]
;   maggies - observed 12, and 22 micron fluxes, as output
;     by WISE_TO_MAGGIES [2,NGAL]
;   ivarmaggies - corresponding inverse variances [4,NGAL]
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   chary - use the Chary & Elbaz (2001) models (default)
;   dale - use the Dale & Helou (2002) models
;   rieke - use the Rieke et al. (2009) models
;
; OUTPUTS: 
;   lir - infrared luminosity for each object [L_sun, NGAL]
;
; OPTIONAL OUTPUTS:
;   err_lir - *statistical* uncertainty on LIR [L_sun, NGAL]
;   model_indx - fractional index number of the best-fitting model 
;     [NGAL] 

; COMMENTS:
;   By default the code uses the Chary & Elbaz (2001) models, unless
;   either /RIEKE or /DALE are set.  The statistical uncertainty on
;   ERR_LIR are derived using a simple Monte Carlo method.  This
;   routine should be called using various models to estimate the
;   *systematic* errors.  
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Apr 19, UCSD
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

function wise2lir_model_maggies, model, zref, dlum=dlum
; compute model photometry on a fiducial redshift grid    
    nmodel = n_elements(model.lir)
    nzref = n_elements(zref)
; k_projection_table can't deal with the units so we have to
; loop and call k_project_filters    
;   k_projection_table, rmatrix, model.flux, model.wave, $
;     zref, (wise_filterlist())[2:3]
;   model_maggies = dblarr(4,nmodel,nzref)
    model_maggies = dblarr(2,nmodel,nzref)
    for ii = 0, nmodel-1 do begin
       for jj = 0, nzref-1 do begin
          model_maggies[*,ii,jj] = k_project_filters(k_lambda_to_edges($
            model.wave*(1.0+zref[jj])),model.flux[*,ii]/(4.0*!dpi*dlum[jj]^2)/$
            (1.0+zref[jj]),filterlist=(wise_filterlist())[2:3])
       endfor
    endfor
return, model_maggies
end

function im_wise2lir, redshift, maggies, ivarmaggies, chi2=chi2, $
  model_indx=model_indx, err_lir=err_lir, chary=chary, $
  dale=dale, rieke=rieke, debug=debug, zlog=zlog

    common com_wise2lir, chary_maggies, dale_maggies, $
      rieke_maggies
    
    ngal = n_elements(redshift)
    if (ngal eq 0L) then begin
       doc_library, 'im_wise2lir'
       return, -1
    endif

;   if (ngal ne n_elements(f24)) then begin
;      splog, 'Dimensions of REDSHIFT and F24 do not agree'
;      return, -1
;   endif
;   if (n_elements(f24_err) ne 0L) then begin
;      if (ngal ne n_elements(f24_err)) then begin
;         splog, 'Dimensions of REDSHIFT and F24_ERR do not agree'
;         return, -1
;      endif
;   endif

; define the fiducial redshift grid
    zrefmin = min(redshift)
    zrefmax = max(redshift)
    if zrefmin eq zrefmax then begin
       zrefmin = zrefmin*0.95
       zrefmax = zrefmax*1.05
    endif
    nzz = (ceil((zrefmax-zrefmin)/0.01)>50)<200
    zref = range(zrefmin,zrefmax,nzz,log=zlog)
    dlum = dluminosity(zref,/cm)
    
; read the relevant set of models and build the fiducial table of
; photometry; use Chary & Elbaz (2001) by default
    if (keyword_set(chary) eq 0) and (keyword_set(dale) eq 0) and $
      (keyword_set(rieke) eq 0) then chary = 1
    
; Chary & Elbaz (2001)    
    if keyword_set(chary) then begin
       model = read_01chary()
       if (n_elements(chary_maggies) eq 0) then $
         chary_maggies = wise2lir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = chary_maggies
    endif
; Dale & Helou (2002)
    if keyword_set(dale) then begin
       model = read_02dale()
       if (n_elements(dale_maggies) eq 0) then $
         dale_maggies = wise2lir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = dale_maggies
    endif
; Rieke et al. (2009)
    if keyword_set(rieke) then begin
       model = read_09rieke()
       if (n_elements(rieke_maggies) eq 0) then $
         rieke_maggies = wise2lir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = rieke_maggies
    endif
    nmodel = n_elements(model.lir)
    npix = n_elements(model.wave)

    filt = (wise_filterlist())[2:3]
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)
    
; interpolate the models at the redshifts of interest
    zindx = findex(zref,redshift)
    model_maggies = interpolate(model_maggies_grid,zindx)
    idlum = interpolate(dlum,zindx)
    
; give the 12- and 22-micron bands equal statistical weight
    new_ivarmaggies = ivarmaggies
    new_ivarmaggies[0,*] = max(ivarmaggies,dim=1)
    new_ivarmaggies[1,*] = max(ivarmaggies,dim=1)
    
;; give the 12- and 22-micron bands equal statistical weight; ignore
;; the 3.3- and 4.6-micron fluxes
;    new_ivarmaggies = ivarmaggies
;    new_ivarmaggies[0:1,*] = 0.0
;;   new_ivarmaggies[2,*] = 0.0
;    new_ivarmaggies[2,*] = max(ivarmaggies,dim=1)
;    new_ivarmaggies[3,*] = max(ivarmaggies,dim=1)
    
; need to loop, unfortunately...  derive the uncertainties on L(IR)
; and L(24) using a simple Monte Carlo method
    model_indx = fltarr(ngal)-1.0
    lir = model_indx*0.0-1.0
    err_lir = model_indx*0.0-1.0
    chi2 = model_indx*0.0-1.0

    nmonte = 100
    nfine = 20

    lir_monte = fltarr(nmonte)
    for igal = 0L, ngal-1 do begin
       nmaggies = abs(maggies[*,igal]*1.0D) ; need absolute value to deal with negative fluxes correctly
       nivarmaggies = new_ivarmaggies[*,igal]*1.0D
       if (total(nivarmaggies gt 0) ge 1) then begin
;      if (total((nmaggies gt 0) and (nivarmaggies gt 0)) ge 1) then begin
; compute and minimize chi2 (see MFINDCHI2MIN)
          vmodelmaggies = model_maggies[*,*,igal]
          vmaggies = rebin(reform(nmaggies,nfilt,1),nfilt,nmodel)
          vivarmaggies = rebin(reform(nivarmaggies,nfilt,1),nfilt,nmodel)
          vchi2 = total(vivarmaggies*(vmaggies-vmodelmaggies)^2,1,/double)

; add a floor to prevent negative minima, which can occur          
          chi2offset = 100
          vchi2 = vchi2+chi2offset 

          xx = findgen(nmodel)
          minx = min(xx,max=maxx)
          y2 = spl_init(xx,vchi2)
          x2 = findgen(nmodel*nfine+1)/(nmodel*nfine)*(maxx-minx)+minx
          chi2fit = spl_interp(xx,vchi2,y2,x2)

; use a lower-order interpolation scheme if the chi^2 goes negative          
          if (min(chi2fit) lt 0.0) then chi2fit = interpol(vchi2,xx,x2,/quad)
          if (min(chi2fit) lt 0.0) then chi2fit = interpol(vchi2,xx,x2)
          minchi2 = min(chi2fit,fitplace)

          modelindx1 = x2[fitplace] 
          model_indx[igal] = modelindx1 ; fractional index number
          lir[igal] = interpolate(model.lir,modelindx1)

          if (minchi2 le 0) then message, 'Negative chi^2'
          minchi2 = minchi2-chi2offset
          chi2[igal] = minchi2
          
; debugging QAplots
          if keyword_set(debug) then begin
             splog, minchi2, modelindx1, lir[igal]

             djs_plot, x2, chi2fit, psym=-6, sym=0.5, xsty=3, ysty=3, /ylog, $
               xrange=modelindx1+5*[-1,1]
             djs_oplot, xx, vchi2, psym=-6, sym=3, color='red'
             djs_oplot, modelindx1*[1,1], 10^!y.crange, color='dark green'
             djs_oplot, !x.crange, minchi2*[1,1], color='dark green'
             cc = get_kbrd(1)
             
             modelflam = interpolate(model.flux,modelindx1,/grid)/(4.0*!dpi*idlum[igal]^2)/(1.0+redshift[igal])
             good = where(modelflam gt 0.0,ngood)
             modelwave = model.wave[good]*(1.0+redshift[igal])
             modelmab = -2.5*alog10(modelflam[good]*rebin(reform(modelwave,ngood,1),$
               ngood,nmodel)^2/im_light(/ang))-48.6

             djs_plot, weff/1D4, -2.5*alog10(maggies[*,igal]), psym=7, color='yellow', $
               sym=3, xr=[1,500], /xlog, yr=[max(modelmab),min(modelmab)], xsty=3, ysty=3
             djs_oplot, weff/1D4, -2.5*alog10(interpolate(vmodelmaggies,modelindx1)), $
               psym=6, sym=3
;            djs_plot, weff*(1.0+redshift[igal])/1D4, -2.5*alog10(maggies[*,igal]), psym=7, color='yellow', $
;              sym=3, xr=[1,500], /xlog, yr=[max(modelmab),min(modelmab)], xsty=3, ysty=3
;            djs_oplot, weff*(1.0+redshift[igal])/1D4, -2.5*alog10(interpolate(vmodelmaggies,modelindx1)), $
;              psym=6, sym=3
             djs_oplot, modelwave/1D4, modelmab, color='orange'

;            plot, findgen(nmodel), vchi2, xsty=3, ysty=3, psym=6, /ylog, $
;              position=[2.0,max(modelmab)+2,9.0,min(modelmab)-2], /noerase, $
;              xrange=[0.0,nmodel], yrange=[minmax(vchi2)], /data
;            niceprint, nmaggies, nivarmaggies
             cc = get_kbrd(1)
          endif

;      if arg_present(err_lir) then begin
;         maggies_monte = maggies[igal] + errmaggies[igal]*randomn(seed,nmonte)
;         for jj = 0, nmonte-1 do begin
;            model_indx1 = interpol(findgen(nmodel),$ ; fractional index number
;              maggies_monte[jj]/model_maggies[*,igal],1.0,/quad)
;            lir_monte[jj] = interpolate(model.lir,model_indx1)
;         endfor
;         err_lir[igal] = djsig(lir_monte)
       endif 
    endfor

return, lir
end
