pro plot_redshift_vs_sfr_mass_density, postscript=postscript, pdf=pdf
; jm08apr08nyu - written, based on PLOT_04HOPKINS
;                (2004) SFRD vs redshift plot

; default plotting variables - note that POSTSCRIPT over-rides PDF

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    charsize = 2.0
    textcolor1 = 'white'
    pspath = getenv('PAPERSPATH')+'/literature/'
    datapath = pspath+'data/'

    if keyword_set(pdf) then begin
       postscript = 0L
       pspath = pspath+'keynote/' ; for keynote presentations
       textcolor1 = 'white'
       charsize = 2.2
       postthick1 = 8.0
       postthick2 = 6.0
       postthick3 = 10.0
    endif
    if keyword_set(postscript) then begin
       textcolor1 = 'black'
       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 6.0
    endif

    timelabel1 = [1.0,3.0,5.0,7.0,10.0,12.0] ; [Gyr]
    timelabel2 = [1.0,3.0,5.0,7.0]

; ---------------------------------------------------------------------------
; REDSHIFT VS SFR AND STELLAR MASS DENSITY    
; ---------------------------------------------------------------------------

    xsize = 8.5 & ysize = 7.9
    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.6, height=3.0*[1,1], $
      xmargin=[1.5,0.4], ymargin=[0.8,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'redshift_vs_sfr_mass_density'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    zrange = [-0.1,7.3]
    xrange = alog10(1+zrange)
    xtitle = 'log (1 + z)'

    zaxis = findgen((7.0-0.0)/0.01+1)*0.01+0.0
    xaxis = alog10(1+zaxis)
    ageaxis = getage(0.0)-getage(zaxis)
    
; upper panel - SFR density
    
    h = rsex(datapath+'04hopkins.sex')
    uv = where(strmatch(h.indicator,'*UV*',/fold))
    ha = where(strmatch(h.indicator,'*Ha*',/fold) or $
      strmatch(h.indicator,'*Hb*',/fold) or $
      strmatch(h.indicator,'*OII*',/fold))
    ir = where(strmatch(h.indicator,'*IR*',/fold))
    rad = where(strmatch(h.indicator,'*RADIO*',/fold) or $
      strmatch(h.indicator,'*xray*',/fold))

    uvcolor = 'royal blue'
    hacolor = 'forest green'
    ircolor = 'firebrick'
    radcolor = 'orchid'
    
    yrange = [-2.3,-0.2]
;   ytitle = 'log \rho_{SFR} !c (M'+sunsymbol()+' yr^{-1} Mpc^{-3})'
;   ytitle = 'log \rho_{SFR}'
    ytitle = 'log \rho_{SFR} (M'+sunsymbol()+' yr^{-1} Mpc^{-3})'
    
    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle=textoidl(ytitle), xtickname=replicate(' ',10), $
      charsize=charsize, charthick=postthick2, xsty=9, ysty=1, xrange=xrange, $
       position=pos[*,0], yrange=yrange, color=fsc_color(textcolor1,150)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=charsize, $
      charthick=postthick2, xtitle='Lookback Time (Gyr)', color=fsc_color(textcolor1,100), $
      xtickv=alog10(1.0+getredshift(getage(0.0)-timelabel1)), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')

; UV    
    plotsym, 8, 0.8, /fill
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_lo/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_lo, ps=8, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_hi/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_hi, ps=8, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; Ha
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_lo/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_lo, ps=8, color=fsc_color(hacolor,10), errcolor=fsc_color(hacolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_hi/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_hi, ps=8, color=fsc_color(hacolor,10), errcolor=fsc_color(hacolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; IR
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_lo/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_lo, ps=8, color=fsc_color(ircolor,10), errcolor=fsc_color(ircolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_hi/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_hi, ps=8, color=fsc_color(ircolor,10), errcolor=fsc_color(ircolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; Radio
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_lo/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_lo, ps=8, color=fsc_color(radcolor,10), errcolor=fsc_color(radcolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_hi/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_hi, ps=8, color=fsc_color(radcolor,10), errcolor=fsc_color(radcolor,10), $
      thick=postthick1, errthick=postthick1, /hibar

    sfrd = 10.0^h.sfrd
    sfrderr = total([[h.sfrderr_lo],[h.sfrderr_hi]],2)/2.0*h.sfrd*alog(10.0)*0.0+0.1 ; NO WEIGHTING!
    pp = fltarr(4)+1.0
    pp = mpfitexpr('(P[0]+P[1]*X)/(1.0+(X/P[2])^P[3])',h.z,sfrd,sfrderr,pp,/quiet)

;   pp = [0.0166,0.1848,1.9474,2.6316] ; Cole et al. 2001
;   pp = [0.014,0.11,1.4,2.2]          ; Wilkins et al. 2008, for h=0.7
    splog, 'SFR density coefficients ', pp

    sfrd_fit = (pp[0]+pp[1]*zaxis)/(1.0+(zaxis/pp[2])^pp[3])
    oplot, xaxis, alog10(sfrd_fit), line=0, thick=postthick3, color=fsc_color(textcolor1)
    
; lower panel - stellar mass density

    w = rsex(datapath+'08wilkins.sex') & nw = n_elements(w)
    zz = fltarr(nw)
    zzerr = fltarr(nw)
    for ii = 0L, nw-1L do begin
       zz[ii] = mean([w[ii].zmin,w[ii].zmax])
       zzerr[ii] = stddev([w[ii].zmin,w[ii].zmax])
    endfor
       
    yrange = [-4.3,-2.01]
    ytitle = 'log \Omega_{*}'
    
    plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize, charthick=postthick2, xsty=1, ysty=1, xrange=xrange, $
      position=pos[*,1], yrange=yrange, color=fsc_color(textcolor1,150)

    plotsym, 8, 0.8, /fill
    oploterror, alog10(1+zz), alog10(w.omega), zzerr/(1+zz)/alog(10.0), $
      w.omega_err/w.omega/alog(10.0), ps=8, thick=postthick1, errthick=postthick1, $
      color=fsc_color('orange',10), errcolor=fsc_color('orange',10)

    cc = fltarr(3)+1.0
    cc = mpfitexpr('(P[0]*exp(-P[1]*X^P[2]))',zz,w.omega,w.omega_err*0.0+1.0,cc,/quiet)
;   cc = [0.0023,0.68,1.2] ; Wilkins et al. (2008)
;   print, cc
    
    omegastar_fit = cc[0]*exp(-cc[1]*zaxis^cc[2])
    oplot, xaxis, alog10(omegastar_fit), line=5, thick=postthick3, color=fsc_color(textcolor1)

; integrate the SFRD curve and overplot on the mass-density plot

    rho_crit = 1.35D11 ; M_sun/Mpc^3 (h=0.7)
    
    mstar_predict = im_integral(ageaxis,sfrd_fit*1D9,ageaxis,$ ; M_sun/Mpc^3
      replicate(max(ageaxis),n_elements(ageaxis)))
    omegastar_predict = (1.0-0.28)*mstar_predict/rho_crit
    
;   oplot, xaxis, alog10(omegastar_predict), line=2, thick=postthick1, color=fsc_color(textcolor1)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ---------------------------------------------------------------------------
; REDSHIFT VS SFR DENSITY    
; ---------------------------------------------------------------------------

    xsize = 8.5 & ysize = 8.4
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.6,0.4], ymargin=[0.8,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'redshift_vs_sfr_density'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)

    h = rsex(datapath+'04hopkins.sex')
    uv = where(strmatch(h.indicator,'*UV*',/fold))
    ha = where(strmatch(h.indicator,'*Ha*',/fold) or $
      strmatch(h.indicator,'*Hb*',/fold) or $
      strmatch(h.indicator,'*OII*',/fold))
    ir = where(strmatch(h.indicator,'*IR*',/fold))
    rad = where(strmatch(h.indicator,'*RADIO*',/fold) or $
      strmatch(h.indicator,'*xray*',/fold))

    uvcolor = 'royal blue'
    hacolor = 'forest green'
    ircolor = 'firebrick'
    radcolor = 'orange'
    
    xrange = [-0.05,0.92]
    yrange = [-2.3,-0.2]
    xtitle = 'log (1 + z)'
    ytitle = 'log \rho_{SFR} (M'+sunsymbol()+' yr^{-1} Mpc^{-3} h_{70}^{3})'
    
    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize, charthick=postthick2, xsty=9, ysty=1, xrange=xrange, $
      position=pos, yrange=yrange, color=fsc_color(textcolor1,150)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=charsize, $
      charthick=postthick2, xtitle='Lookback Time (Gyr)', color=fsc_color(textcolor1,100), $
      xtickv=alog10(1.0+getredshift(getage(0.0)-timelabel1)), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')

; UV    
    plotsym, 8, 0.8, /fill
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_lo/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_lo, ps=8, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[uv].z), h[uv].sfrd, h[uv].zerr_hi/(1+h[uv].z)/alog(10.0), $
      h[uv].sfrderr_hi, ps=8, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; Ha
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_lo/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_lo, ps=8, color=fsc_color(hacolor,10), errcolor=fsc_color(hacolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[ha].z), h[ha].sfrd, h[ha].zerr_hi/(1+h[ha].z)/alog(10.0), $
      h[ha].sfrderr_hi, ps=8, color=fsc_color(hacolor,10), errcolor=fsc_color(hacolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; IR
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_lo/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_lo, ps=8, color=fsc_color(ircolor,10), errcolor=fsc_color(ircolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[ir].z), h[ir].sfrd, h[ir].zerr_hi/(1+h[ir].z)/alog(10.0), $
      h[ir].sfrderr_hi, ps=8, color=fsc_color(ircolor,10), errcolor=fsc_color(ircolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
; Radio
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_lo/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_lo, ps=8, color=fsc_color(radcolor,10), errcolor=fsc_color(radcolor,10), $
      thick=postthick1, errthick=postthick1, /lobar
    oploterror, alog10(1+h[rad].z), h[rad].sfrd, h[rad].zerr_hi/(1+h[rad].z)/alog(10.0), $
      h[rad].sfrderr_hi, ps=8, color=fsc_color(radcolor,10), errcolor=fsc_color(radcolor,10), $
      thick=postthick1, errthick=postthick1, /hibar
    
;; overplot some pegase models
;
;    peg1 = mrdfits(getenv('PEGASE_HR_SFHGRID_DIR')+'/MEASURE/salp_tau_001.0Gyr.info.fits',1,silent=0)
;    zform = 1.5 & tform = getage(zform)
;    peg_z = getredshift(peg1.age/1D3+tform) ; same time array for all Pegase models
;    djs_oplot, alog10(1+peg_z), peg1.sfr*1E11, line=5, thick=2
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ---------------------------------------------------------------------------
; REDSHIFT VS SFR DENSITY, SPLIT BY MASS
; ---------------------------------------------------------------------------

    xsize = 8.5 & ysize = 8.4
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.6,0.4], ymargin=[0.8,1.1], xpage=xsize, ypage=ysize, $
      position=pos, /normal

    psname = 'redshift_vs_sfr_density_bymass'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xsize, ysize=ysize, silent=keyword_set(pdf)
    
    h = rsex(datapath+'04hopkins.sex')

    xrange = [-0.05,1.1]
    yrange = [-3.5,-0.7]
    xtitle = 'Redshift'
    ytitle = 'log \rho_{SFR} (M'+sunsymbol()+' yr^{-1} Mpc^{-3})'
;   ytitle = 'log \rho_{SFR} (M'+sunsymbol()+' yr^{-1} Mpc^{-3} h_{70}^{3})'
    
    plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize, charthick=postthick2, xsty=9, ysty=1, xrange=xrange, $
      position=pos, yrange=yrange, color=fsc_color(textcolor1,150)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=charsize, $
      charthick=postthick2, xtitle='Lookback Time (Gyr)', color=fsc_color(textcolor1,100), $
      xtickv=getredshift(getage(0.0)-timelabel2), xticks=n_elements(timelabel2)-1L, $
      xtickname=string(timelabel2,format='(I0)')

    zbin1 = [0.3,0.5,0.65,0.9]
    zbin1_err = [0.1,0.1,0.15,0.1]
    zbin2 = [0.01,0.3,0.5,0.65,0.9]
    zbin2_err = [0.0,0.1,0.1,0.15,0.1]
    zbin3 = [0.01,0.3,0.5,0.65,0.9]
    zbin3_err = [0.0,0.1,0.1,0.15,0.1]

    massbin1 = [-2.22,-1.94,-1.72,-1.68]
    massbin1_err_up = [0.13,0.12,0.09,0.09]
    massbin1_err_lo = [0.18,0.17,0.11,0.11]
    
    massbin2 = [-2.26,-1.93,-1.71,-1.49,-1.45]
    massbin2_err_up = [0.08,0.10,0.07,0.05,0.06]
    massbin2_err_lo = [0.10,0.14,0.08,0.05,0.07]

    massbin3 = [-3.01,-2.72,-2.41,-2.25,-2.28]
    massbin3_err_up = [0.13,0.25,0.16,0.12,0.14]
    massbin3_err_lo = [0.19,0.63,0.26,0.17,0.22]

    totalmass = alog10(total(10^[[massbin1],[massbin2[1:4]],[massbin3[1:4]]],2))
    totalmass_err_up = sqrt(total([[massbin1_err_up],[massbin2_err_up[1:4]],$
      [massbin3_err_up[1:4]]]^2.0,2))
    totalmass_err_lo = sqrt(total([[massbin1_err_lo],[massbin2_err_lo[1:4]],$
      [massbin3_err_lo[1:4]]]^2.0,2))
    
    plotsym, 8, 2.5, /fill ; 9<logM<10
    oploterror, zbin1, massbin1, zbin1_err, massbin1_err_up, /hibar, $
      ps=-8, color=fsc_color('blue',8), errcolor=fsc_color('blue',8), $
      thick=postthick1, errthick=postthick1, line=3
    oploterror, zbin1, massbin1, zbin1_err, massbin1_err_lo, /lobar, $
      ps=-8, color=fsc_color('blue',8), errcolor=fsc_color('blue',8), $
      thick=postthick1, errthick=postthick1, line=3

    plotsym, 5, 2.5, /fill ; 10<logM<11
    oploterror, zbin2, massbin2, zbin2_err, massbin2_err_up, /hibar, $
      ps=-8, color=fsc_color('forest green',10), errcolor=fsc_color('forest green',10), $
      thick=postthick1, errthick=postthick1, line=1
    plotsym, 5, 2.5, /fill ; 10<logM<11
    oploterror, zbin2, massbin2, zbin2_err, massbin2_err_lo, /lobar, $
      ps=-8, color=fsc_color('forest green',10), errcolor=fsc_color('forest green',10), $
      thick=postthick1, errthick=postthick1, line=1

    plotsym, 0, 2.5, /fill ; logM>11
    oploterror, zbin3, massbin3, zbin3_err, massbin3_err_up, /hibar, $
      ps=-8, color=fsc_color('red',9), errcolor=fsc_color('red',9), $
      thick=postthick1, errthick=postthick1, line=5
    plotsym, 0, 2.5, /fill ; logM>11
    oploterror, zbin3, massbin3, zbin3_err, massbin3_err_lo, /lobar, $
      ps=-8, color=fsc_color('red',9), errcolor=fsc_color('red',9), $
      thick=postthick1, errthick=postthick1, line=5

; Hopkins (2004)    
    
;   plotsym, 8, 0.3, /fill
;   oploterror, h.z, h.sfrd, h.zerr, h.sfrderr, ps=8, $
;     color=fsc_color('orange',11), errcolor=fsc_color('orange',11), $
;     thick=1.0, errthick=1.0

    plotsym, 8, 2.5, /fill
    oploterror, zbin1, totalmass, zbin1_err, totalmass_err_up, /hibar, $
      ps=-8, color=fsc_color(textcolor1,8), errcolor=fsc_color(textcolor1,8), $
      thick=postthick1, errthick=postthick1
    oploterror, zbin1, totalmass, zbin1_err, totalmass_err_lo, /lobar, $
      ps=-8, color=fsc_color(textcolor1,8), errcolor=fsc_color(textcolor1,8), $
      thick=postthick1, errthick=postthick1

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop
    
return
end
