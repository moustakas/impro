pro age_vs_redshift_polynomial, write=write
; jm06feb23uofa - derive the polynomial that relates age given the
;                 redshift 

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7 ; cosmological parameters

    zmax = 7.0 & zmin = 0.0 & dz = 0.05
    zobj = findgen((zmax-zmin)/dz+1)*dz+zmin
;   print, zobj
    age = getage(zobj)

    coeff = poly_fit(zobj,age,6,yfit=yfit,/double)
    print, coeff
    plot, zobj, age, ps=4, syms=2
    oplot, zobj, yfit, thick=2
;   niceprint, age, yfit

    stats = im_stats(100.0*(age-yfit)/age,/verbose)

return
end
    
