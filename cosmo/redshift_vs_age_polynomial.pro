pro redshift_vs_age_polynomial, write=write
; jm06apr19uofa - derive the polynomial that relates redshift given
;                 the age

    red, h100=0.7, omega_0=0.3, omega_lambda=0.7 ; cosmological parameters

    agemax = getage(0.0) & agemin = 0.0 & dage = 0.01
    age = findgen((agemax-agemin)/dage+1)*dage+agemin
    nage = n_elements(age)
    print, age
    zobj = fltarr(nage)
    for i = 0L, nage-1L do zobj[i] = getredshift(age[i])

    coeff = poly_fit(age,zobj,6,yfit=yfit,/double)
    print, coeff
    plot, age, zobj, ps=4, syms=2
    oplot, age, yfit, thick=2
;   niceprint, age, zobj, yfit

    stats = im_stats(100.0*(zobj-yfit)/zobj,/verbose)

return
end
    
