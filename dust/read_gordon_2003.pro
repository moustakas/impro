function read_gordon_2003, wave, smc=smc, lmc2=lmc2, avglmc=avglmc
; jm04jan26uofa - written
; jm09aug20ucsd - force k(lambda) to be positive at long wavelengths 

    if (not keyword_set(smc)) and (not keyword_set(lmc2)) and (not keyword_set(avglmc)) then begin
       print, 'Syntax - k_lambda = read_gordon_2003(wave,/smc,/lmc2,/avglmc)'
       return, -1L
    endif
    
    path = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='dust')
    file = '2003_gordon.dat'

    if (file_test(path+file) eq 0L) then begin
       splog, 'Gordon et al. (2003) data file not found.'
       return, -1L
    endif

    readcol, path+file, l, x, Al_AV_smc, err, Al_AV_lmc2, err, $
      Al_Av_avglmc, err, comment='#', /silent

    R_V_smc    = 2.74
    R_V_lmc2   = 2.76
    R_V_avglmc = 3.41

    if keyword_set(smc) then begin
       R_V = R_V_smc
       good = where(Al_AV_smc ne -9.999,ngood)
       Al_AV = Al_AV_smc[good]
       l = l[good]
    endif
    if keyword_set(lmc2) then begin
       R_V = R_V_lmc2
       good = where(Al_AV_lmc2 ne -9.999,ngood)
       Al_AV = Al_AV_lmc2[good]
       l = l[good]
    endif
    if keyword_set(avglmc) then begin
       R_V = R_V_avglmc
       good = where(Al_AV_avglmc ne -9.999,ngood)
       Al_AV = Al_AV_avglmc[good]
       l = l[good]
    endif
    
    lambda = 1D4*l
    k_lambda = R_V*Al_AV

    if (n_elements(wave) ne 0L) then k_lambda = interpol(k_lambda,lambda,wave)>0.0 else wave = lambda
    
return, k_lambda
end    
