function init_isedfit, ngal, nfilt, sfhgrid, sfhgrid_paramfile=sfhgrid_paramfile, $
  isedfit_post=isedfit_post
; ISEDFIT support routine - initialize the output structure 

    params = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ndraw = isedfit_ndraw() ; number of random draws

; compute the maximum number of bursts, if any
    if (params.pburst le 0.0) then nmaxburst = 0 else $
      nmaxburst = ceil((params.maxage-params.minage)/params.pburstinterval)

    burstarray1 = -1.0
    burstarray2 = -1.0
    if (nmaxburst gt 1) then begin
       burstarray1 = fltarr(nmaxburst)-1.0
       burstarray2 = dblarr(nmaxburst)-1D
    endif

    isedfit1 = {$
      isedfit_id:      -1L,$    ; unique ID number
      zobj:           -1.0,$    ; redshift
      maggies:     fltarr(nfilt),$ ; observed maggies
      ivarmaggies: fltarr(nfilt),$ ; corresponding inverse variance
      bestmaggies: fltarr(nfilt)}  ; best-fitting model photometry

; best-fit values (at the chi2 minimum); see BUILD_ISEDFIT_SFHGRID
    best = {$
      imf:              '',$
      chunkindx:        -1,$
      modelindx:       -1L,$
      ageindx:          -1,$

      tau:            -1.0,$
      Z:              -1.0,$
      ebv:            -1.0,$
      mu:             -1.0,$
      nburst:            0,$
      tauburst:        -1D,$ ; burst truncation time scale
      tburst:    burstarray2,$
      dtburst:   burstarray2,$
      fburst:    burstarray1,$

      mass:           -1.0,$ 
      age:            -1.0,$
      sfr:            -1.0,$ ; instantaneous
      sfr100:         -1.0,$ ; 100 Myr
      b100:           -1.0,$ ; averaged over the previous 100 Myr
      chi2:            1E6}  ; chi2 minimum

; median quantities and PDF quantiles
    qmed = {$
      mass_avg:     -1.0,$
      age_avg:      -1.0,$
      sfr_avg:      -1.0,$ ; instantaneous
      sfr100_avg:   -1.0,$ ; 100 Myr
      b100_avg:     -1.0,$
      tau_avg:      -1.0,$
      Z_avg:        -1.0,$
      ebv_avg:      -1.0,$
      mu_avg:       -1.0,$

      mass_err:     -1.0,$
      age_err:      -1.0,$
      sfr_err:      -1.0,$
      sfr100_err:   -1.0,$
      b100_err:     -1.0,$
      tau_err:      -1.0,$
      Z_err:        -1.0,$
      ebv_err:      -1.0,$
      mu_err:       -1.0,$

      mass_50:     -1.0,$
      age_50:      -1.0,$
      sfr_50:      -1.0,$ ; instantaneous
      sfr100_50:   -1.0,$ ; 100 Myr
      b100_50:     -1.0,$
      tau_50:      -1.0,$
      Z_50:        -1.0,$
      ebv_50:      -1.0,$
      mu_50:       -1.0,$

      mass_eff_err:   -1.0,$
      age_eff_err:    -1.0,$
      sfr_eff_err:    -1.0,$
      sfr100_eff_err: -1.0,$
      b100_eff_err:   -1.0,$
      tau_eff_err:    -1.0,$
      Z_eff_err:      -1.0,$
      ebv_eff_err:    -1.0,$
      mu_eff_err:     -1.0}

    isedfit = struct_addtags(struct_addtags(best,qmed),isedfit1)
    isedfit = replicate(temporary(isedfit),ngal)

; initialize the posterior distribution structure
    isedfit_post = {$
      draws:  lonarr(ndraw)-1,$
      mass:   fltarr(ndraw)-1}
;     Z:      fltarr(ndraw),$
;     age:    fltarr(ndraw),$
;     tau:    fltarr(ndraw),$
;     ebv:    fltarr(ndraw),$
;     mu:     fltarr(ndraw),$
;     sfr:    fltarr(ndraw),$
;     sfr100: fltarr(ndraw),$
;     b100:   fltarr(ndraw)}
    isedfit_post = replicate(temporary(isedfit_post),ngal)
    
return, isedfit
end
