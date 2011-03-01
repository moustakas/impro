pro fit_smc_lmc

; SMC bar

    readcol, 'smcbar_ext.dat', x, smc, smcerr, skipline=2, /silent
    lambda = 1D4/x  ; Angstrom
    npts = n_elements(lambda)
    
; fit to three different regions

    uv = where(lambda lt 3500.0)
    opt = where((lambda ge 3500.0) and (lambda lt 1E4))
    nir = where(lambda ge 1E4)
    
    uvcoeff = poly_fit(lambda[uv],smc[uv],4,measure_errors=smcerr[uv],$
      yfit=uvfit,chisq=uvchisq)
    optcoeff = poly_fit(lambda[opt],smc[opt],2,measure_errors=smcerr[opt],$
      yfit=optfit,chisq=optchisq)
    nircoeff = poly_fit(lambda[nir],smc[nir],2,measure_errors=smcerr[nir],$
      yfit=nirfit,chisq=nirchisq)

    ploterror, lambda, smc, smcerr, ps=3, xsty=3, ysty=3
    oplot, lambda[uv], uvfit, thick=2.0
    oplot, lambda[opt], optfit, thick=2.0
    oplot, lambda[nir], nirfit, thick=2.0

    coeff = poly_fit(lambda,smc,6,measure_errors=smcerr,$
      yfit=fit,chisq=chisq)

    bkpt = [500.0,1000.0,2000.0,3000.0,3500.0,4000.0,5000.0,6000.0,7000.0,$
      8000.0,9000.0,10000.0,15000.0,20000.0,25000.0]
    sset = bspline_iterfit(lambda,smc,bkpt=bkpt,nord=3,yfit=yfit)
    
return
end    
