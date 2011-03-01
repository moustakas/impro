pro test_im_convolve_sfh
; jm10jan20ucsd - test my code
    
; read the SSP and test tau-model
    minwave = 3500.0 & maxwave = 9000.0
    ssp = im_read_bc03(minwave=minwave,maxwave=maxwave,bc03_extras=essp)

;; --------------------------------------------------
;; test - add one burst
;    infosfh = {tau: 5.0, nburst: 0}
;    f1 = im_convolve_sfh(ssp,infosfh=infosfh,sfh=sfh1,time=time1)
;
;;   infosfh = {tau: 5.0, nburst: 1, tburst: 4D, dtauburst: 0.25D, fburst: 0.4}
;;   f2 = im_convolve_sfh(ssp,infosfh=infosfh,sfh=sfh2,time=time2)
;
;    infosfh = {tau: 5.0, nburst: 1, tburst: 4D, dtauburst: 0.25D, fburst: 0.4}
;    time2 = build_isedfit_agegrid(infosfh,nage=100,minage=min(ssp.age)/1D9,maxage=max(ssp.age)/1D9)
;    sfh2 = 1.22*isedfit_reconstruct_sfh(infosfh,age=time2)
;;   sfh2 = interpol(sfh2,time2,time1)
;;   f2 = im_convolve_sfh(ssp,sfh=sfh2,time=time1)
;    f2 = im_convolve_sfh(ssp,sfh=sfh2,time=time2)
;
;    f21 = interpolate(f2,findex(time2,time1))
;    
;    djs_plot, ssp.wave, f1[*,5]
;    djs_oplot, ssp.wave, f21[*,5], color='red'
;    
;    
;    djs_plot, time1, sfh1, psym=-6
;    djs_oplot, time2, sfh2, color='green', psym=-6
;
;;   infosfh = {tau: 5.0, nburst: 2, tburst: [3.5,5], dtauburst: [0.25,0.4], fburst: [0.5,0.2]}
    
; --------------------------------------------------
; test - compare against BC03 and make a QAplot
    infosfh = {tau: 5.0, nburst: 0}
    time = ssp.age/1D9
    sfh = isedfit_reconstruct_sfh(infosfh,age=time)
    outflux = im_convolve_sfh(ssp,sfh=sfh,time=time,$
      mstar=essp.m_,cspmstar=cspmstar)

    tau = im_read_bc03(isedpath='~/',isedfile='tau5test.ised',$
      minwave=minwave,maxwave=maxwave,bc03_extras=etau)
    tauflux = interpolate(tau.flux,findex(tau.age,time*1D9),/grid)

    niceprint, cspmstar, interpolate(etau.m_,findex(tau.age,time*1D9))

    psfile = '~/test_convolve_sfh.ps'
    im_plotconfig, 6, pos, psfile=psfile
    for ii = 0, n_elements(time)-1 do begin
       djs_plot, ssp.wave, tauflux[*,ii], xsty=3, ysty=3, position=pos[*,0], $
         xtickname=replicate(' ',10), xtitle='', ytitle='Flux'
       djs_oplot, ssp.wave, outflux[*,ii], color='red'
       legend, 'Age = '+strtrim(time[ii],2)+' Gyr', /right, /top, box=0
       djs_plot, ssp.wave, (tauflux[*,ii]-outflux[*,ii])/outflux[*,ii], $
         position=pos[*,1], /noerase, xsty=3, ysty=3, yrange=[-0.1,0.1]
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    
    
return
end
    
