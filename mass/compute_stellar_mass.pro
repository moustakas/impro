;+
; NAME:
;       COMPUTE_STELLAR_MASS()
;
; PURPOSE:
;       Compute stellar mass-to-light ratios and stellar masses using
;       the Bell et al. color-based method. 
;
; CALLING SEQUENCE:
;       stellar_mass = compute_stellar_mass(photo)
;
; INPUTS:
;       photo - input photometry
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       stellar_mass - mass-to-light ratios and stellar masses for a
;                      variety of color and luminosity combinations
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       READ_01BELL()
;
; INTERNAL SUPPORT ROUTINES:
;       COMPUTE_MASS(), COMPUTE_ML()
;
; COMMENTS:
;       We adopt the Salpeter (1955) IMF.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 22, U of A
;       jm04aug01uofa - add 0.15 dex to convert from a "diet"
;                       Salpeter IMF to the Salpeter (1955) IMF
;       jm05jan06uofa - separately compute the mass-to-light ratios;
;                       if synthesized rest-frame magnitudes exist
;                       (eg., from spectra), then also compute the
;                       synthesized mass-to-light ratios
;       jm05may17uofa - completely streamlined and optimized
;-

function compute_mass, photo, ml, ml_err, lumband

    ngalaxy = n_elements(photo)
    tags = tag_names(photo)

    lumgood = tag_exist(photo,lumband,index=lumindx)
    lumerrgood = tag_exist(photo,lumband+'_err',index=lumerrindx)
    
    result = {$
      mass:       -999.0, $
      mass_err:   -999.0}
    result = replicate(result,ngalaxy)

    if (lumindx[0] ne -1L) then begin

       good = where((ml gt -900) and (photo.(lumindx) gt -900.0) and $
         (photo.(lumerrindx) gt -900.0),ngood)

       if (ngood ne 0L) then begin

          lum = 10.0^photo[good].(lumindx)
          lum_err = alog(10.0)*lum*photo[good].(lumerrindx)

          mass = lum*ml[good]
          mass_err = sqrt((lum*ml_err[good])^2 + (lum_err*ml[good])^2 )
             
          result[good].mass = alog10(mass)
          result[good].mass_err = mass_err/mass/alog(10.0)

       endif
          
    endif
          
return, result
end

function compute_ml, photo, band1, band2, coeff

    imf_dietsalpeter_to_salpeter = +0.15

    ngalaxy = n_elements(photo)
    tags = tag_names(photo)
    
    b1good = tag_exist(photo,band1,index=b1indx)
    b2good = tag_exist(photo,band2,index=b2indx)

    b1errgood = tag_exist(photo,band1+'_err',index=b1errindx)
    b2errgood = tag_exist(photo,band2+'_err',index=b2errindx)

    result = {$
      ml:       -999.0, $
      ml_err:   -999.0}
    result = replicate(result,ngalaxy)

    if (b1indx[0] ne -1L) and (b2indx[0] ne -1L) and (b1errindx[0] ne -1L) and (b2errindx[0] ne -1L) then begin

       good = where((photo.(b1indx) gt -900.0) and (photo.(b2indx) gt -900.0) and $
         (photo.(b1errindx) gt -900.0) and (photo.(b2errindx) gt -900.0),ngood)

       if (ngood ne 0L) then begin

          color = photo[good].(b1indx)-photo[good].(b2indx)
          color_err = sqrt(photo[good].(b1errindx)^2+photo[good].(b2errindx)^2)
          
          result[good].ml = 10.0^(coeff[0] + coeff[1]*color + imf_dietsalpeter_to_salpeter) ; note IMF offset
          result[good].ml_err = result[good].ml*coeff[1]*color_err*alog(10.0)

       endif
       
    endif else begin

       splog, 'Insufficient information to compute a mass-to-light '+$
         'ratio using bands: '+strupcase(band1)+', '+strupcase(band2)+'.'
       
    endelse
    
return, result
end

function compute_stellar_mass, photo

    ngalaxy = n_elements(photo)
    if (ngalaxy eq 0L) then begin
       print, 'Syntax - mass = compute_stellar_mass(photo)'
       return, -1L
    endif

; read the Bell et al. coefficients

    bell01 = read_01bell(bell03=bell03)
    
; initialize the output structure    

    stellar_mass = {$
;     ml_bv_b:             -999.0, $
;     ml_bv_b_err:         -999.0, $
;     ml_bv_b_2003:        -999.0, $ ; using the Bell et al. (2003) coefficients
;     ml_bv_b_2003_err:    -999.0, $
      ml_bv_v:             -999.0, $
      ml_bv_v_err:         -999.0, $
;     ml_br_b:             -999.0, $
;     ml_br_b_err:         -999.0, $
      ml_br_r:             -999.0, $
      ml_br_r_err:         -999.0, $
;     ml_br_b_2003:        -999.0, $
;     ml_br_b_2003_err:    -999.0, $
;     ml_br_r_2003:        -999.0, $
;     ml_br_r_2003_err:    -999.0, $
;     ml_vh_v:             -999.0, $
;     ml_vh_v_err:         -999.0, $
;     ml_vh_h:             -999.0, $
;     ml_vh_h_err:         -999.0, $
;     ml_vk_v:             -999.0, $
;     ml_vk_v_err:         -999.0, $
;     ml_vk_k:             -999.0, $
;     ml_vk_k_err:         -999.0, $
;     ml_sdss_gr_g:        -999.0, $
;     ml_sdss_gr_g_err:    -999.0, $
      ml_sdss_gr_r:        -999.0, $
      ml_sdss_gr_r_err:    -999.0, $
;     ml_sdss_ri_r:        -999.0, $
;     ml_sdss_ri_r_err:    -999.0, $
;     ml_sdss_ri_i:        -999.0, $
;     ml_sdss_ri_i_err:    -999.0, $

;     mass_bv_b:           -999.0, $
;     mass_bv_b_err:       -999.0, $
;     mass_bv_b_2003:      -999.0, $ ; using the Bell et al. (2003) coefficients
;     mass_bv_b_2003_err:  -999.0, $
      mass_bv_v:           -999.0, $
      mass_bv_v_err:       -999.0, $
;     mass_br_b:           -999.0, $
;     mass_br_b_err:       -999.0, $
      mass_br_r:           -999.0, $
      mass_br_r_err:       -999.0, $
;     mass_br_b_2003:      -999.0, $
;     mass_br_b_2003_err:  -999.0, $
;     mass_br_r_2003:      -999.0, $
;     mass_br_r_2003_err:  -999.0, $
;     mass_vh_v:           -999.0, $
;     mass_vh_v_err:       -999.0, $
;     mass_vh_h:           -999.0, $
;     mass_vh_h_err:       -999.0, $
;     mass_vk_v:           -999.0, $
;     mass_vk_v_err:       -999.0, $
;     mass_vk_k:           -999.0, $
;     mass_vk_k_err:       -999.0, $
;     mass_sdss_gr_g:      -999.0, $
;     mass_sdss_gr_g_err:  -999.0, $
      mass_sdss_gr_r:      -999.0, $
      mass_sdss_gr_r_err:  -999.0};, $
;     mass_sdss_ri_r:      -999.0, $
;     mass_sdss_ri_r_err:  -999.0, $
;     mass_sdss_ri_i:      -999.0, $
;     mass_sdss_ri_i_err:  -999.0}

    stellar_mass = replicate(stellar_mass,ngalaxy)

; ---------------------------------------------------------------------------    
; (B-V), B - 2001 coefficients
; ---------------------------------------------------------------------------    

;   coeff = [bell01[0].ab,bell01[0].bb] ; Bell & de Jong (2001)
;   
;   ml = compute_ml(photo,'B','V',coeff)
;   stellar_mass.ml_bv_b     = ml.ml
;   stellar_mass.ml_bv_b_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'B_lum')
;   stellar_mass.mass_bv_b     = mass.mass
;   stellar_mass.mass_bv_b_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-V), B - 2003 coefficients
; ---------------------------------------------------------------------------    

;   coeff = [bell03[9].ag,bell03[9].bg] ; NOTE: ag=ab, bg=bg, etc.
;   
;   ml = compute_ml(photo,'B','V',coeff)
;   stellar_mass.ml_bv_b_2003     = ml.ml
;   stellar_mass.ml_bv_b_2003_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'B_lum')
;   stellar_mass.mass_bv_b_2003     = mass.mass
;   stellar_mass.mass_bv_b_2003_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-V), V
; ---------------------------------------------------------------------------    

    coeff = [bell01[0].av,bell01[0].bv]
    
    ml = compute_ml(photo,'B','V',coeff)
    stellar_mass.ml_bv_v     = ml.ml
    stellar_mass.ml_bv_v_err = ml.ml_err

    mass = compute_mass(photo,ml.ml,ml.ml_err,'V_lum')
    stellar_mass.mass_bv_v     = mass.mass
    stellar_mass.mass_bv_v_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-R), B - 2001 coefficients
; ---------------------------------------------------------------------------    

;   coeff = [bell01[1].ab,bell01[1].bb]
;   
;   ml = compute_ml(photo,'B','R',coeff)
;   stellar_mass.ml_br_b     = ml.ml
;   stellar_mass.ml_br_b_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'B_lum')
;   stellar_mass.mass_br_b     = mass.mass
;   stellar_mass.mass_br_b_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-R), R - 2001 coefficients
; ---------------------------------------------------------------------------    

    coeff = [bell01[1].ar,bell01[1].br]
    
    ml = compute_ml(photo,'B','R',coeff)
    stellar_mass.ml_br_r     = ml.ml
    stellar_mass.ml_br_r_err = ml.ml_err

    mass = compute_mass(photo,ml.ml,ml.ml_err,'R_lum')
    stellar_mass.mass_br_r     = mass.mass
    stellar_mass.mass_br_r_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-R), B - 2003 coefficients
; ---------------------------------------------------------------------------    

;   coeff = [bell03[10].ag,bell03[10].bg]
;   
;   ml = compute_ml(photo,'B','R',coeff)
;   stellar_mass.ml_br_b_2003     = ml.ml
;   stellar_mass.ml_br_b_2003_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'B_lum')
;   stellar_mass.mass_br_b_2003     = mass.mass
;   stellar_mass.mass_br_b_2003_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (B-R), R - 2003 coefficients
; ---------------------------------------------------------------------------    

;   coeff = [bell03[10].ai,bell03[10].bi] ; NOTE: ai=ar, bi=br, etc.
;   
;   ml = compute_ml(photo,'B','R',coeff)
;   stellar_mass.ml_br_r_2003     = ml.ml
;   stellar_mass.ml_br_r_2003_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'R_lum')
;   stellar_mass.mass_br_r_2003     = mass.mass
;   stellar_mass.mass_br_r_2003_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (V-H), V
; ---------------------------------------------------------------------------    

;   coeff = [bell01[4].av,bell01[4].bv]
;   
;   ml = compute_ml(photo,'V','H',coeff)
;   stellar_mass.ml_vh_v     = ml.ml
;   stellar_mass.ml_vh_v_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'V_lum')
;   stellar_mass.mass_vh_v     = mass.mass
;   stellar_mass.mass_vh_v_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (V-H), H
; ---------------------------------------------------------------------------    

;   coeff = [bell01[4].ah,bell01[4].bh]
;   
;   ml = compute_ml(photo,'V','H',coeff)
;   stellar_mass.ml_vh_h     = ml.ml
;   stellar_mass.ml_vh_h_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'H_lum')
;   stellar_mass.mass_vh_h     = mass.mass
;   stellar_mass.mass_vh_h_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (V-K), V
; ---------------------------------------------------------------------------    

;   coeff = [bell01[5].av,bell01[5].bv]
;   
;   ml = compute_ml(photo,'V','K',coeff)
;   stellar_mass.ml_vk_v     = ml.ml
;   stellar_mass.ml_vk_v_err = ml.ml_err
;
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'V_lum')
;   stellar_mass.mass_vk_v     = mass.mass
;   stellar_mass.mass_vk_v_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (V-K), K
; ---------------------------------------------------------------------------    

;   coeff = [bell01[5].ak,bell01[5].bk]
;   
;   ml = compute_ml(photo,'V','K',coeff)
;   stellar_mass.ml_vk_k     = ml.ml
;   stellar_mass.ml_vk_k_err = ml.ml_err
;   
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'K_lum')
;   stellar_mass.mass_vk_k     = mass.mass
;   stellar_mass.mass_vk_k_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (g-r), g
; ---------------------------------------------------------------------------    

;   coeff = [bell03[4].ag,bell03[4].bg]
;   
;   ml = compute_ml(photo,'sdss_g','sdss_r',coeff)
;   stellar_mass.ml_sdss_gr_g     = ml.ml
;   stellar_mass.ml_sdss_gr_g_err = ml.ml_err
;   
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'sdss_g_lum')
;   stellar_mass.mass_sdss_gr_g     = mass.mass
;   stellar_mass.mass_sdss_gr_g_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (g-r), r
; ---------------------------------------------------------------------------    

    coeff = [bell03[4].ar,bell03[4].br]
    
    ml = compute_ml(photo,'sdss_g','sdss_r',coeff)
    stellar_mass.ml_sdss_gr_r     = ml.ml
    stellar_mass.ml_sdss_gr_r_err = ml.ml_err
    
    mass = compute_mass(photo,ml.ml,ml.ml_err,'sdss_r_lum')
    stellar_mass.mass_sdss_gr_r     = mass.mass
    stellar_mass.mass_sdss_gr_r_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (r-i), r
; ---------------------------------------------------------------------------    

;   coeff = [bell03[7].ar,bell03[7].br]
;   
;   ml = compute_ml(photo,'sdss_r','sdss_i',coeff)
;   stellar_mass.ml_sdss_ri_r     = ml.ml
;   stellar_mass.ml_sdss_ri_r_err = ml.ml_err
;   
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'sdss_r_lum')
;   stellar_mass.mass_sdss_ri_r     = mass.mass
;   stellar_mass.mass_sdss_ri_r_err = mass.mass_err

; ---------------------------------------------------------------------------    
; (r-i), i
; ---------------------------------------------------------------------------    

;   coeff = [bell03[7].ai,bell03[7].bi]
;   
;   ml = compute_ml(photo,'sdss_r','sdss_i',coeff)
;   stellar_mass.ml_sdss_ri_i     = ml.ml
;   stellar_mass.ml_sdss_ri_i_err = ml.ml_err
;   
;   mass = compute_mass(photo,ml.ml,ml.ml_err,'sdss_i_lum')
;   stellar_mass.mass_sdss_ri_i     = mass.mass
;   stellar_mass.mass_sdss_ri_i_err = mass.mass_err

return, stellar_mass
end    
