pro test_simple_kcorrect

    ff = 'sdss_'+['u0','g0','r0','i0','z0']+'.par'
;   ff = 'bessell_'+['B','R','I']+'.par'
    nfilt = n_elements(ff)

    band_shift = 0.1
    
    dist = 10.0*3.085678D18 ; 10 pc fiducial distance [cm]

    sed = im_read_bc03(age=[0.1,1.0,5.0,10.0])
    z = [0.25,0.33,0.5,0.95]
    mass = [0.1,1,5,10.0]*1D10
    ngal = n_elements(z)
    npix = n_elements(sed.wave)
    
    obsmaggies = dblarr(nfilt,ngal)
    obsmaggies_ivar = dblarr(nfilt,ngal)

    sedwave = dblarr(npix,ngal)
    sedflux = dblarr(npix,ngal)
    for jj = 0L, ngal-1L do begin
       sedwave[*,jj] = sed.wave*(1+z[jj])
       sedflux[*,jj] = mass[jj]*3.826D33*sed.flux[*,jj]/(1+z[jj])/$
         (4.0*!pi*dluminosity(z[jj],/cm)^2.0)
       
       obsmaggies[*,jj] = k_project_filters(k_lambda_to_edges(reform(sedwave[*,jj])),$
         reform(sedflux[*,jj]),filterlist=ff,/silent)
       obsmaggies_err = obsmaggies[*,jj]/20.0D
       obsmaggies[*,jj] = obsmaggies[*,jj] + randomn(seed,nfilt)*obsmaggies_err
       obsmaggies_ivar[*,jj] = 1.0/obsmaggies_err^2.0
    endfor

    restwave = sed.wave
    restflux = 3.826D33*sed.flux/(4*!dpi*dist^2.0) ; [erg/s/cm2/A/M_sun]

    kcorr1 = im_simple_kcorrect(z,obsmaggies,obsmaggies_ivar,$
      ff,restwave,restflux,absmag=absmag1,rmaggies=rmaggies1,$
      chi2=chi21,scale=scale1,/sdss,/silent,band_shift=band_shift)
    mass1 = scale1*(dluminosity(z,/cm)/dist)^2.0 ; stellar mass

    kcorr2 = im_kcorrect(z,obsmaggies,obsmaggies_ivar,ff,mass=mass2,$
      absmag=absmag2,rmaggies=rmaggies2,chi2=chi2,coeffs=coeffs,/sdss,$
      band_shift=band_shift)

    niceprint, mass1, mass2, mass

    col = ['red','green','blue','cyan','purple','orange']
    plot, [0], [0], xrange=[-17,-23], yrange=[-17,-23], xsty=3, ysty=3
    oplot, !x.crange, !y.crange, line=0
    for ii = 0L, nfilt-1L do djs_oplot, absmag1[ii,*], absmag2[ii,*], ps=4, $
      color=col[ii]
    
stop    
    
return
end
    
