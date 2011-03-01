;+
; NAME:
;   IM_F24TOLIR()
;
; PURPOSE:
;   Given a redshift and an observed 24-micron flux, compute the total
;   infrared luminosity, L(IR)=L(8-1100), based on various infrared
;   SED models.
;
; INPUTS: 
;   redshift - redshift for each object [NGAL]
;   f24 - observed 24-micron flux [NGAL, mJy]
;
; OPTIONAL INPUTS: 
;   f24_err - 1-sigma uncertainty on F24; mandatory if ERR_LIR or
;     ERR_L24 are requested [NGAL, mJy] 
;   f24_limit - limiting flux [scalar, mJy]; if provided, then ZMAX is
;     computed 
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
;   err_lir - statistical uncertainty on LIR [L_sun, NGAL]
;   err_l24 - statistical uncertainty on L24 [L_sun, NGAL]
;   model_indx - fractional index number of the best-fitting model 
;     [NGAL] 
;   zmax - if F24_LIMIT is provided, then ZMAX is computed as the
;     redshift at which the 24-micron flux equals F24_LIMIT [NGAL]

; COMMENTS:
;   By default the code uses the Chary & Elbaz (2001) models, unless
;   either /RIEKE or /DALE are set.  The statistical uncertainty on
;   ERR_LIR and/or ERR_L24 are derived using a simple Monte Carlo
;   method.  This routine should be called using various models to
;   estimate the *systematic* errors. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 17, UCSD 
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

function f24tolir_model_maggies, model, zref, dlum=dlum
; compute model photometry on a fiducial redshift grid    
    nmodel = n_elements(model.lir)
    nzref = n_elements(zref)
    model_maggies = dblarr(nmodel,nzref)
    for ii = 0, nmodel-1 do begin
       for jj = 0, nzref-1 do begin
          model_maggies[ii,jj] = k_project_filters(k_lambda_to_edges($
            model.wave*(1.0+zref[jj])),model.flux[*,ii]/(4.0*!dpi*dlum[jj]^2)/$
            (1.0+zref[jj]),filterlist='spitzer_mips_24.par')
       endfor
    endfor
return, model_maggies
end

function im_f24tolir, redshift, f24, f24_err=f24_err, $
  model_indx=model_indx, err_lir=err_lir, l24=l24, $
  err_l24=err_l24, f24_limit=f24_limit, zmax=zmax, $
  chary=chary, dale=dale, rieke=rieke

    common com_f24tolir, chary_maggies, dale_maggies, $
      rieke_maggies
    
    ngal = n_elements(redshift)
    if (ngal eq 0L) then begin
       doc_library, 'im_f24tolir'
       return, -1
    endif
    if (ngal ne n_elements(f24)) then begin
       splog, 'Dimensions of REDSHIFT and F24 do not agree'
       return, -1
    endif
    if (n_elements(f24_err) ne 0L) then begin
       if (ngal ne n_elements(f24_err)) then begin
          splog, 'Dimensions of REDSHIFT and F24_ERR do not agree'
          return, -1
       endif
    endif
    if arg_present(err_lir) and (n_elements(f24_err) eq 0L) then begin
       splog, 'Computing ERR_LIR requires F24_ERR'
       return, -1
    endif

; convert the input photometry to maggies
    factor = 1D-3*1D-23*10.0^(0.4*48.6) ; [mJy] --> [maggies]
    maggies = f24*factor
    if (n_elements(f24_err) ne 0) then errmaggies = f24_err*factor
    if (n_elements(f24_limit) ne 0) then maggies_limit = f24_limit*factor
    
; define the fiducial redshift grid    
    zrefmin = 0.02
    zrefmax = 7.0
    zref = im_array(zrefmin,zrefmax,0.02)
    dlum = dluminosity(zref,/cm)
    
; read the relevant set of models and build the fiducial table of
; photometry; use Chary & Elbaz (2001) by default
    if (keyword_set(chary) eq 0) and (keyword_set(dale) eq 0) and $
      (keyword_set(rieke) eq 0) then chary = 1
    
; Chary & Elbaz (2001)    
    if keyword_set(chary) then begin
       model = read_01chary()
       if (n_elements(chary_maggies) eq 0) then $
         chary_maggies = f24tolir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = chary_maggies
    endif
; Dale & Helou (2002)
    if keyword_set(dale) then begin
       model = read_02dale()
       if (n_elements(dale_maggies) eq 0) then $
         dale_maggies = f24tolir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = dale_maggies
    endif
; Rieke et al. (2009)
    if keyword_set(rieke) then begin
       model = read_09rieke()
       if (n_elements(rieke_maggies) eq 0) then $
         rieke_maggies = f24tolir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = rieke_maggies
    endif
    nmodel = n_elements(model.lir)

; interpolate the models at the redshifts of interest    
    model_maggies = interpolate(model_maggies_grid,$
      findex(zref,redshift),/grid)
    
; need to loop, unfortunately...  derive the uncertainties on L(IR)
; and L(24) using a simple Monte Carlo method
    model_indx = fltarr(ngal)
    lir = model_indx*0.0
    l24 = model_indx*0.0
    err_lir = model_indx*0.0
    err_l24 = model_indx*0.0
    if (n_elements(f24_limit) ne 0) then zmax = fltarr(ngal)
    
    nmonte = 100
    lir_monte = fltarr(nmonte)
    l24_monte = fltarr(nmonte)
    for ii = 0L, ngal-1L do begin
       model_indx[ii] = interpol(findgen(nmodel),$ ; fractional index number
         maggies[ii]/model_maggies[*,ii],1.0,/quad)
       lir[ii] = interpolate(model.lir,model_indx[ii])
       l24[ii] = interpolate(model.l24,model_indx[ii])
       if arg_present(err_lir) or arg_present(err_l24) then begin
          maggies_monte = maggies[ii] + errmaggies[ii]*randomn(seed,nmonte)
          for jj = 0, nmonte-1 do begin
             model_indx1 = interpol(findgen(nmodel),$ ; fractional index number
               maggies_monte[jj]/model_maggies[*,ii],1.0,/quad)
             lir_monte[jj] = interpolate(model.lir,model_indx1)
             l24_monte[jj] = interpolate(model.l24,model_indx1)
          endfor
          err_lir[ii] = djsig(lir_monte)
          err_l24[ii] = djsig(l24_monte)
       endif
; compute ZMAX for this object *assuming the same template*!
       if (n_elements(f24_limit) ne 0) then zmax[ii] = $
         interpol(zref,interpolate(transpose(model_maggies_grid),model_indx[ii])/$
         maggies_limit,1.0,/quad)
    endfor

return, lir
end
