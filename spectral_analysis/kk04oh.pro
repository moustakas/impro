function kk04oh, r23, r23_err, o32, o32_err
; jm06apr18uofa - compute the oxygen abundance based on the KK04 R23
;                 calibration; excised from IM_ABUNDANCE() for testing
;                 with SIMULATE_R23OH

    nmonte = 100L
    light = 2.99792458D10 ; speed of light [cm/s]
    
    nobj = n_elements(r23)
    if (nobj eq 0L) then begin
       print, 'Syntax - '
       return, -1L
    endif
    
    result = {zstrong_o32:                     0.0, $
              zstrong_o32_err:                 0.0, $
              zstrong_r23:                     0.0, $
              zstrong_r23_err:                 0.0, $
              zstrong_12oh_kk04_r23_upper:     0.0, $
              zstrong_12oh_kk04_r23_upper_err: 0.0, $
              zstrong_12oh_kk04_logu_upper:         0.0, $
              zstrong_12oh_kk04_r23_lower:     0.0, $
              zstrong_12oh_kk04_r23_lower_err: 0.0, $
              zstrong_12oh_kk04_logu_lower: 0.0}
    result = replicate(result,nobj)

    result.zstrong_r23     = r23
    result.zstrong_r23_err = r23_err
    result.zstrong_o32     = o32
    result.zstrong_o32_err = o32_err

    x = r23 & xerr = r23_err
    y = o32 & yerr = o32_err

    ymonte = rebin(reform(y,1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*rebin(reform(yerr,1,nobj),nmonte,nobj)
    xmonte = rebin(reform(x,1,nobj),nmonte,nobj)+randomn(seed,nmonte,nobj)*rebin(reform(xerr,1,nobj),nmonte,nobj)

; upper branch       

    logoh_upper = 9.0

    maxiter = 10L
    for j = 0L, maxiter do begin ; iterate
       logq = (32.81D - 1.153*y^2 + logoh_upper*(-3.396 - 0.025*y + 0.1444*y^2)) / $
         (4.603D - 0.3119*y - 0.163*y^2 + logoh_upper*(-0.48 + 0.0271*y + 0.02037*y^2))
       logoh_upper = 9.72D - 0.777*x - 0.951*x^2 - 0.072*x^3 - 0.811*x^4 - $
         logq*(0.0737 - 0.0713*x - 0.141*x^2 + 0.0373*x^3 - 0.058*x^4)
    endfor

    result.zstrong_12oh_kk04_r23_upper = temporary(logoh_upper)
    result.zstrong_12oh_kk04_logu_upper = logq - alog10(light)

    logqmonte = rebin(reform(logq,1,nobj),nmonte,nobj)
    logoh_upper_monte = 9.72D - 0.777*xmonte - 0.951*xmonte^2 - 0.072*xmonte^3 - 0.811*xmonte^4 - $
      logqmonte*(0.0737 - 0.0713*xmonte - 0.141*xmonte^2 + 0.0373*xmonte^3 - 0.058*xmonte^4)
    for ii = 0L, nobj-1L do result[ii].zstrong_12oh_kk04_r23_upper_err = stddev(logoh_upper_monte[*,ii])

; lower branch       
    
    logoh_lower = 7.9
    
    for j = 0L, maxiter do begin ; iterate

       logq = (32.81D - 1.153*y^2 + logoh_lower*(-3.396 - 0.025*y + 0.1444*y^2)) / $
         (4.603D - 0.3119*y - 0.163*y^2 + logoh_lower*(-0.48 + 0.0271*y + 0.02037*y^2))
       logoh_lower = 9.40D + 4.65*x - 3.17*x^2 - logq*(0.272 + 0.547*x - 0.513*x^2)

    endfor

    result.zstrong_12oh_kk04_r23_lower = temporary(logoh_lower)
    result.zstrong_12oh_kk04_logu_lower = logq - alog10(light)
    
    logqmonte = rebin(reform(logq,1,nobj),nmonte,nobj)
    logoh_lower_monte = 9.40D + 4.65*xmonte - 3.17*xmonte^2 - logqmonte*(0.272 + 0.547*xmonte - 0.513*xmonte^2)
    for ii = 0L, nobj-1L do result[ii].zstrong_12oh_kk04_r23_lower_err = stddev(logoh_lower_monte[*,ii])
       
return, result
end
