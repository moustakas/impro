pro benchmark_isedfit

    nmodel = long(1E5)

    nfilt = 10
    maggies = make_array(nfilt,/nozero)
    ivarmaggies = make_array(nfilt,/nozero)

    modelmaggies = make_array(nfilt,nmodel,/nozero)

    vmaggies = rebin(reform(maggies,nfilt,1),nfilt,nmodel)
    vivarmaggies = rebin(reform(ivarmaggies,nfilt,1),nfilt,nmodel)
    vmodelmaggies = reform(modelmaggies,nfilt,nmodel)

    t0 = systime(1)
    vmass = total(reform((ivarmaggies*maggies),1,nfilt)#vmodelmaggies,1)/$
      total(reform(ivarmaggies,1,nfilt)#vmodelmaggies^2.0,1)
    vchi2 = total(vivarmaggies*(vmaggies-rebin(reform(vmass,1,nmodel),$
      nfilt,nmodel)*vmodelmaggies)^2.0,1)
    print, nmodel, systime(1)-t0
    
return
end
    
