pro build_all_isedfit_sfhgrids
; jm09may22nyu - build all the various SFH grids; only call this
;   routine if you know what you're doing!

; Chabrier IMF
    t0 = systime(1)
    build_isedfit_sfhgrid, synthmodels='bc03', imf='chab', $
      sfhgrid=[1,2,3], redcurve=[0,1,2,3], /clobber
    splog, 'Total time = ', (systime(1)-t0)/60.0

stop
;; Salpeter IMF
;    t0 = systime(1)
;    build_isedfit_sfhgrid, synthmodels='bc03', imf='salp', $
;      sfhgrid=[1,2,3], redcurve=[0,1,2,3], clobber=0
;    splog, 'Total time = '+(systime(1)-t0)/60.0
    
return
end
    
