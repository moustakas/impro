function pfunction, x, p, R23=R23

    a1 = p[0]
    a2 = p[1]
    a3 = p[2]
    b1 = p[3]
    b2 = p[4]
    b3 = p[5]
    b0 = p[6]

    Phard = x
    
    model = (R23 - a1 - a2*Phard - a3*Phard^2) / (b1 + b2*Phard + b3*Phard^2 - b0*R23)

return, model
end

pro minimize_pcoefficients
; jm04jul08uofa

    hii = read_hii_regions(/nolog)
    indx = where((hii.oii_h_beta gt -900.0) and (hii.oiii_h_beta gt -900.0) and $
      (hii.ZT_log12OH gt -900.0) and (hii.ZT_log12oh gt 8.2),nindx)
    hii = hii[indx]

    R2 = hii.oii_h_beta
    R3 = hii.oiii_h_beta
    R23 = R2 + R3

    P = R3 / R23

    OH12 = hii.ZT_log12OH
    OH12_err = hii.ZT_log12OH_err

    lo = where(OH12 lt 8.2,comp=hi)
    djs_plot, P[lo], R3[lo], ps=4, xsty=3, ysty=3, color='green', $
      yrange=[0,14]
    djs_oplot, P[hi], R3[hi], ps=7, color='red'

    djs_plot, P, OH12, ps=4, xsty=3, ysty=3
    djs_plot, R23, OH12, ps=4, xsty=3, ysty=3

    z = OH12 # (fltarr(nindx)+1.0)
    x = P # replicate(1,nindx)
    y = replicate(1,nindx) # R23
    
    acoeff = djs_sfit(z,x,y,3,4,yfit=yfit)
    
    
    
    
    lo = where(P lt 0.5,comp=hi)
    djs_plot, alog10(R23[lo]), OH12[lo], ps=4, xsty=3, ysty=3, color='green', $
      xrange=[0.0,1.2], yrange=[7.5,9.5]
    djs_oplot, alog10(R23[hi]), OH12[hi], ps=7, color='red'

    nparams = 7L
    parinfo = {$
      parname:     'Coeff', $
      value:          0.0D, $
      fixed:            0L, $
      limited:     [0L,0L], $
      limits:  [0.0D,0.0D], $
      step:              0, $ ; automatic
      relstep:           0, $ ; automatic
      mpside:            0, $ ; automatic
      mpmaxstep:         0, $ ; automatic
      tied:             '', $
      mpprint:           1}
    parinfo = replicate(parinfo,nparams)
    parinfo.value = [-54.2,-59.45,-7.31,+6.07,+6.71,+0.371,-0.243]

    functargs = {R23: R23}

    coeff = mpfitfun('pfunction',P,OH12,OH12_err,/autoderivative,$
      bestnorm=bestnorm,covar=covar,dof=dof,ftol=ftol,yfit=yfit,$
      functargs=functargs,gtol=gtol,maxiter=maxiter,nfev=nfev,$
      nfree=nfree,niter=niter,npegged=npegged,parinfo=parinfo,$
      perror=perror,/quiet,resdamp=0,status=status,xtol=xtol)

    mpchi2 = sqrt(bestnorm/dof)
    splog, 'MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
    if (status EQ 5) then splog, 'Warning: Maximum number of iterations reached: ', niter

stop    
    

return
end
    
