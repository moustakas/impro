function read_smc, wave, r_v=r_v
; jm02sep27uofa
; SMC bar

    if n_elements(r_v) eq 0L then r_v = 3.1

    path = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='dust')
    readcol, path+'smcbar_ext.dat', x, Al_AV, skipline=2, /silent

    lambda = 1D4/x  ; Angstrom
    smc = r_v*Al_AV

    if n_elements(wave) ne 0L then begin

      smc = interpol(smc,lambda,wave,/spline) 
      
    endif else wave = lambda
    
return, smc
end    
