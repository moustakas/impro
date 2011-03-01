pro im_panelplot, theta, magloss, label, postscript=postscript

    pagemaker, nx=2, ny=2, position=position, /normal, yspace=0.0, $
      xspace=0.0, xmargin=[1.3,1.0], ymargin=[0.5,1.0]

    yrange = [0.01,max(magloss)]
    
    djs_plot, theta, magloss[*,0,0], xsty=3, ysty=3, xthick=2.0, ythick=2.0, $
      charsize=1.5, charthick=2.0, xrange=xrange, yrange=yrange, $
      position=position[*,0], line=0, thick=3.0, color='blue', $
      xtickname=replicate(' ',10), ytitle='Light Loss (mag)', /ylog
    djs_oplot, theta, magloss[*,0,1], line=1, thick=3.0, color='green'
    djs_oplot, theta, magloss[*,0,2], line=2, thick=3.0, color='red'
    legend, label[*,0], /left, /top, box=0, charsize=1.5, charthick=2.0

    djs_plot, theta, magloss[*,1,0], xsty=3, ysty=11, xthick=2.0, ythick=2.0, $
      charsize=1.5, charthick=2.0, xrange=xrange, yrange=yrange, $
      position=position[*,1], line=0, thick=3.0, color='blue', /noerase, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), /ylog
    axis, /yaxis, yrange=(10^(0.4*yrange)-1.0)*100.0, charsize=1.5, charthick=2.0, $
      ythick=2.0, ytitle='Light Loss (%)', ysty=3
    djs_oplot, theta, magloss[*,1,1], line=1, thick=3.0, color='green'
    djs_oplot, theta, magloss[*,1,2], line=2, thick=3.0, color='red'
    legend, label[*,1], /left, /top, box=0, charsize=1.5, charthick=2.0

    djs_plot, theta, magloss[*,2,0], xsty=3, ysty=3, xthick=2.0, ythick=2.0, $
      charsize=1.5, charthick=2.0, xrange=xrange, yrange=yrange, $
      position=position[*,2], line=0, thick=3.0, color='blue', /noerase, $
      xtitle='\Delta\theta (degrees)', ytitle='Light Loss (mag)', /ylog
    djs_oplot, theta, magloss[*,2,1], line=1, thick=3.0, color='green'
    djs_oplot, theta, magloss[*,2,2], line=2, thick=3.0, color='red'
    legend, label[*,2], /left, /top, box=0, charsize=1.5, charthick=2.0

    djs_plot, theta, magloss[*,3,0], xsty=3, ysty=11, xthick=2.0, ythick=2.0, $
      charsize=1.5, charthick=2.0, xrange=xrange, yrange=yrange, $
      position=position[*,3], line=0, thick=3.0, color='blue', /noerase, $
      ytickname=replicate(' ',10), xtitle='\Delta\theta (degrees)', /ylog
    axis, /yaxis, yrange=(10^(0.4*yrange)-1.0)*100.0, charsize=1.5, charthick=2.0, $
      ythick=2.0, ytitle='Light Loss (%)', ysty=3
    djs_oplot, theta, magloss[*,3,1], line=1, thick=3.0, color='green'
    djs_oplot, theta, magloss[*,3,2], line=2, thick=3.0, color='red'
    legend, label[*,3], /left, /top, box=0, charsize=1.5, charthick=2.0
    
    if not keyword_set(postscript) then cc = get_kbrd(1)

return
end

function calculate_seeing, wave, seeing_V
; compute the wavelength-dependent seeing by choosing the seeing at 
; 5000 Angstrom (V-band)
    
    r0 = wave[1]*206265.0/seeing_V   ; [Angstrom]
    r0_vec = r0*(wave/5000.0)^(6./5) ; [Angstrom]

    seeing = 206265*wave/r0_vec ; [arcsec]
    
return, seeing
end    

pro im_lightloss, postscript=postscript
; jm03apr21uofa

    wave = [3727.0,5000.0,6563.0] ; [Angstrom]
    nwave = n_elements(wave)

; --------------------------------------------------
; SLIT WIDTH = 4.5"
; --------------------------------------------------
    
    slit_width = 4.5
    airmass = [1.0,1.2,1.5,2.0]
    nairmass = n_elements(airmass)

    theta = findgen((90.0-0.0)/5.0+1.0)*5.0+0.0 ; [0-90] degrees
    ntheta = n_elements(theta)
    xrange = minmax(theta)

; --------------------------------------------------    
; PAGE 1    
; --------------------------------------------------    

    seeing_V = 1.5 
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 1.5"','4.5" Slit'], $
      ['Airmass = 1.2','Seeing = 1.5"',''], $
      ['Airmass = 1.5','Seeing = 1.5"',''], $
      ['Airmass = 2.0','Seeing = 1.5"',''] ]
      
    loss = fltarr(ntheta,nairmass,nwave)
    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    if keyword_set(postscript) then dfpsplot, 'im_lightloss_4.5.ps', $
      /square, /color, /isolatin1       

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 1
    
; --------------------------------------------------    
; PAGE 2
; --------------------------------------------------    

    seeing_V = 2.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 2.0"','4.5" Slit'], $
      ['Airmass = 1.2','Seeing = 2.0"',''], $
      ['Airmass = 1.5','Seeing = 2.0"',''], $
      ['Airmass = 2.0','Seeing = 2.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2
    
; --------------------------------------------------    
; PAGE 3
; --------------------------------------------------    

    seeing_V = 3.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 3.0"','4.5" Slit'], $
      ['Airmass = 1.2','Seeing = 3.0"',''], $
      ['Airmass = 1.5','Seeing = 3.0"',''], $
      ['Airmass = 2.0','Seeing = 3.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

; --------------------------------------------------    
; PAGE 4
; --------------------------------------------------    

    seeing_V = 4.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 4.0"','4.5" Slit'], $
      ['Airmass = 1.2','Seeing = 4.0"',''], $
      ['Airmass = 1.5','Seeing = 4.0"',''], $
      ['Airmass = 2.0','Seeing = 4.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

    if keyword_set(postscript) then dfpsclose else cc = get_kbrd(1)

; --------------------------------------------------
; SLIT WIDTH = 2.5"
; --------------------------------------------------
    
    slit_width = 2.5
    airmass = [1.0,1.2,1.5,2.0]
    nairmass = n_elements(airmass)

    theta = findgen((90.0-0.0)/5.0+1.0)*5.0+0.0 ; [0-90] degrees
    ntheta = n_elements(theta)
    xrange = minmax(theta)

; --------------------------------------------------    
; PAGE 1    
; --------------------------------------------------    

    seeing_V = 1.5 
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 1.5"','2.5" Slit'], $
      ['Airmass = 1.2','Seeing = 1.5"',''], $
      ['Airmass = 1.5','Seeing = 1.5"',''], $
      ['Airmass = 2.0','Seeing = 1.5"',''] ]
      
    loss = fltarr(ntheta,nairmass,nwave)
    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    if keyword_set(postscript) then dfpsplot, 'im_lightloss_2.5.ps', $
      /square, /color, /isolatin1       

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 1
    
; --------------------------------------------------    
; PAGE 2
; --------------------------------------------------    

    seeing_V = 2.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 2.0"','2.5" Slit'], $
      ['Airmass = 1.2','Seeing = 2.0"',''], $
      ['Airmass = 1.5','Seeing = 2.0"',''], $
      ['Airmass = 2.0','Seeing = 2.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2
    
; --------------------------------------------------    
; PAGE 3
; --------------------------------------------------    

    seeing_V = 3.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 3.0"','2.5" Slit'], $
      ['Airmass = 1.2','Seeing = 3.0"',''], $
      ['Airmass = 1.5','Seeing = 3.0"',''], $
      ['Airmass = 2.0','Seeing = 3.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

; --------------------------------------------------    
; PAGE 4
; --------------------------------------------------    

    seeing_V = 4.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 4.0"','2.5" Slit'], $
      ['Airmass = 1.2','Seeing = 4.0"',''], $
      ['Airmass = 1.5','Seeing = 4.0"',''], $
      ['Airmass = 2.0','Seeing = 4.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

    if keyword_set(postscript) then dfpsclose

; --------------------------------------------------
; SLIT WIDTH = 1.0"
; --------------------------------------------------
    
    slit_width = 1.0
    airmass = [1.0,1.2,1.5,2.0]
    nairmass = n_elements(airmass)

    theta = findgen((90.0-0.0)/5.0+1.0)*5.0+0.0 ; [0-90] degrees
    ntheta = n_elements(theta)
    xrange = minmax(theta)

; --------------------------------------------------    
; PAGE 1    
; --------------------------------------------------    

    seeing_V = 1.5 
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 1.5"','1.0" Slit'], $
      ['Airmass = 1.2','Seeing = 1.5"',''], $
      ['Airmass = 1.5','Seeing = 1.5"',''], $
      ['Airmass = 2.0','Seeing = 1.5"',''] ]
      
    loss = fltarr(ntheta,nairmass,nwave)
    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    if keyword_set(postscript) then dfpsplot, 'im_lightloss_1.0.ps', $
      /square, /color, /isolatin1       

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 1
    
; --------------------------------------------------    
; PAGE 2
; --------------------------------------------------    

    seeing_V = 2.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 2.0"','1.0" Slit'], $
      ['Airmass = 1.2','Seeing = 2.0"',''], $
      ['Airmass = 1.5','Seeing = 2.0"',''], $
      ['Airmass = 2.0','Seeing = 2.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2
    
; --------------------------------------------------    
; PAGE 3
; --------------------------------------------------    

    seeing_V = 3.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 3.0"','1.0" Slit'], $
      ['Airmass = 1.2','Seeing = 3.0"',''], $
      ['Airmass = 1.5','Seeing = 3.0"',''], $
      ['Airmass = 2.0','Seeing = 3.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

; --------------------------------------------------    
; PAGE 4
; --------------------------------------------------    

    seeing_V = 4.0
    seeing = calculate_seeing(wave,seeing_V) ; [arcsec]
    
    label = [ $
      ['Airmass = 1.0','Seeing = 4.0"','1.0" Slit'], $
      ['Airmass = 1.2','Seeing = 4.0"',''], $
      ['Airmass = 1.5','Seeing = 4.0"',''], $
      ['Airmass = 2.0','Seeing = 4.0"',''] ]

    for i = 0L, nwave-1L do for j = 0L, nairmass-1L do for k = 0L, ntheta-1L do $
      loss[k,j,i] = starlightloss(seeing[i],slit_width,theta=theta[k],$
      airmass=airmass[j],wave0=wave[i],max_arcsec=max_arcsec,/noplot,/silent)
    magloss = 2.5*alog10(1.0+loss)

    im_panelplot, theta, magloss, label, postscript=postscript ; PAGE 2

    if keyword_set(postscript) then dfpsclose

return
end    
