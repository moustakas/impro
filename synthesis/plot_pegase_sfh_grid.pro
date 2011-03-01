;+
; NAME:
;       PLOT_PEGASE_SFH_GRID
;
; PURPOSE:
;       Plot quantities of interest in the Pegase-HR SFH grid. 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       A postscript plot is written to
;       getenv('PEGASE_HR_SFHGRID_DIR').   
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       Run after WRITE_PEGASE_SFH_GRID!
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 19, U of A: based on
;         WRITE_SFH_BASE_MODELS 
;       jm07nov13nyu - major updates/rewrite
;
; Copyright (C) 2006-2007, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro plot_pegase_sfh_grid, zform=zform, infall=infall, postscript=postscript

; constants    
    
    basepath = getenv('PEGASE_HR_SFHGRID_DIR')+'/' ; SFH grid path
    gridpath = basepath+'MEASURE/' 

    red, h100=0.7, omega0=0.3, omegal=0.7
    
    if (n_elements(zform) eq 0L) then zform = 5.0

    mgalaxy = 1D11 ; M_sun
    tform = getage(zform)   ; age of the universe when the galaxy formed 
    tuniverse = getage(0.0) ; age of the universe today

    z1 = 0.0       & t1 = getage(z1)
    z2 = 5.0<zform & t2 = getage(z2)
;   z2 = zform & t2 = getage(z2)
    
    lookbackt = findgen((tuniverse-tform)/0.1D)*0.1D ; [0.0-TFORM]
    lookbackz = getredshift(tuniverse-lookbackt)>0.0 ; [0.0-ZFORM]
    nlook = n_elements(lookbackt)

    labelage = t1-[0.0,5.0,10.0,12.0]
    labelz = getredshift(labelage)>0.0
    keep = where(labelz lt zform)
    labelage = labelage[keep]
    labelz = labelz[keep]
    niceprint, labelage, labelz
    
; tau-models with no infall 

;   tau = ['001.0','003.0','008.0','100.0']
;   ntau = n_elements(tau)
;
;   for itau = 0L, ntau-1L do begin
;
;      peg_tau1 = mrdfits(gridpath+'salp_tau_'+tau[itau]+'Gyr.info.fits',1,/silent)
;      pegtime = peg_tau1.age ; [Gyr]
;
;      if (itau eq 0L) then peg_tau = im_empty_structure(peg_tau1,ncopies=nlook)
;
;      for itag = 0L, n_tags(peg_tau1)-1L do if (size(peg_tau1.(itag),/type) eq 7L) then $
;        peg_tau.(itag) = peg_tau1[0].(itag) else peg_tau.(itag) = interpol(peg_tau1.(itag),pegtime,lookbackt*1D3)
;
;      if (n_elements(allpeg_tau) eq 0L) then allpeg_tau = peg_tau else allpeg_tau = [ [allpeg_tau], [peg_tau]]
;      
;   endfor

    if keyword_set(infall) then begin
       
; Kennicutt-Schmidt models with infall 

       peg1 = mrdfits(gridpath+'salp_kennlaw_0.10_infall_001.0Gyr.info.fits',1,/silent)
       peg2 = mrdfits(gridpath+'salp_kennlaw_0.10_infall_003.0Gyr.info.fits',1,/silent)
       peg3 = mrdfits(gridpath+'salp_kennlaw_0.10_infall_005.0Gyr.info.fits',1,/silent)
       peg4 = mrdfits(gridpath+'salp_kennlaw_0.10_infall_008.0Gyr.info.fits',1,/silent)
       peg5 = mrdfits(gridpath+'salp_kennlaw_0.10_infall_100.0Gyr.info.fits',1,/silent)

       tau_suffix = '_{inf}'
       psname = 'plot_salp_kennlaw_infall_zform_'+string(zform,format='(F3.1)')+'.ps'
       
    endif else begin

; tau-models with no infall 

       peg1 = mrdfits(gridpath+'salp_tau_001.0Gyr.info.fits',1,/silent)
       peg2 = mrdfits(gridpath+'salp_tau_003.0Gyr.info.fits',1,/silent)
       peg3 = mrdfits(gridpath+'salp_tau_005.0Gyr.info.fits',1,/silent)
       peg4 = mrdfits(gridpath+'salp_tau_008.0Gyr.info.fits',1,/silent)
       peg5 = mrdfits(gridpath+'salp_tau_999.0Gyr.info.fits',1,/silent)

       tau_suffix = '_{sfr}'
       psname = 'plot_salp_tau_zform_'+string(zform,format='(F3.1)')+'.ps'
       
    endelse
       
    taulines = [0,3,2] & taucolors = ['dodger blue','red','orange']
;   taulines = [0,3,5] & taucolors = ['dodger blue','saddle brown','orange']

    get_element, peg1.age/1D3, [t2,t1]-tform, age_indx
    peg1 = peg1[age_indx[0]:age_indx[1]]
    peg2 = peg2[age_indx[0]:age_indx[1]]
    peg3 = peg3[age_indx[0]:age_indx[1]]
    peg4 = peg4[age_indx[0]:age_indx[1]]
    peg5 = peg5[age_indx[0]:age_indx[1]]
    
;; models with infall 
;    
;    ipeg1 = mrdfits(gridpath+'salp_tau_01.0Gyr_infall_01.0Gyr.fits',1,/silent) ; tau=tauinfall
;    ipeg3 = mrdfits(gridpath+'salp_tau_04.0Gyr_infall_04.0Gyr.fits',1,/silent) ; tau=tauinfall
;    ipeg6 = mrdfits(gridpath+'salp_tau_08.0Gyr_infall_08.0Gyr.fits',1,/silent) ; tau=tauinfall
;
;    ipeg2 = mrdfits(gridpath+'salp_tau_04.0Gyr_infall_01.0Gyr.fits',1,/silent) ; tau>tinfall
;    ipeg5 = mrdfits(gridpath+'salp_tau_08.0Gyr_infall_04.0Gyr.fits',1,/silent) ; tau>tinfall
;    ipeg9 = mrdfits(gridpath+'salp_tau_15.0Gyr_infall_08.0Gyr.fits',1,/silent) ; tau>tinfall
;                                                                              
;    ipeg4 = mrdfits(gridpath+'salp_tau_08.0Gyr_infall_01.0Gyr.fits',1,/silent) ; tau>>tinfall
;    ipeg7 = mrdfits(gridpath+'salp_tau_15.0Gyr_infall_01.0Gyr.fits',1,/silent) ; tau>>tinfall
;    ipeg8 = mrdfits(gridpath+'salp_tau_15.0Gyr_infall_04.0Gyr.fits',1,/silent) ; tau>>tinfall

; plotting variables    
    
    if keyword_set(postscript) then begin
       dfpsplot, basepath+psname, /color, /square
       postthick1 = 5.0
       postthick2 = 3.0
       postthick3 = 4.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick1 = 2.0
       postthick2 = 2.0
       postthick3 = 2.0
    endelse

    plotsym, 0, 1.0, /fill

    mgrange = [-16.5,-25]
    massrange = [7.5,12.5]
    ohrange = [6.8,9.6]
    yeffrange = [-1.8,-2.1]
    grrange = [-0.4,1.3]
    ugrange = [-0.5,2.0]
    
; color-color diagram

    xtitle = '^{0.1}(g - r)'
    ytitle = '^{0.1}(u - g)'

    xrange = grrange
    yrange = ugrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = peg.ugriz[1]-peg.ugriz[2] & y = peg.ugriz[0]-peg.ugriz[1] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = peg.ugriz[1]-peg.ugriz[2] & y = peg.ugriz[0]-peg.ugriz[1] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = peg.ugriz[1]-peg.ugriz[2] & y = peg.ugriz[0]-peg.ugriz[1] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /left, /top, $
      box=0, charsize=2.0, charthick=postthick2
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /right, /bottom, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; luminosity-color diagram

    xtitle = 'M_{0.1g}'
    ytitle = '^{0.1}(g - r)'

    xrange = mgrange
    yrange = grrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /left, /top, $
      box=0, charsize=2.0, charthick=postthick2
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /right, /top, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; luminosity-metallicity diagram

    xtitle = 'M_{0.1g}'
    ytitle = '12 + log (O/H)'

    xrange = mgrange
    yrange = ohrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /right, /top, $
      box=0, charsize=2.0, charthick=postthick2    
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /left, /top, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; luminosity-yeff diagram

    xtitle = 'M_{0.1g}'
    ytitle = 'log (y_{eff})'

    xrange = mgrange
    yrange = yeffrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = alog10(peg.yeff) & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = alog10(peg.yeff) & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = peg.ugriz[1]-2.5*alog10(mgalaxy) & y = alog10(peg.yeff) & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /left, /bottom, $
      box=0, charsize=2.0, charthick=postthick2    
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /right, /bottom, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; mass-color diagram

    xtitle = 'log (M / M_{\odot})'
    ytitle = '^{0.1}(g - r)'

    xrange = massrange
    yrange = grrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = alog10(peg.mstar*mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = alog10(peg.mstar*mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = alog10(peg.mstar*mgalaxy) & y = peg.ugriz[1]-peg.ugriz[2] & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /right, /top, $
      box=0, charsize=2.0, charthick=postthick2    
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /left, /top, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; mass-metallicity diagram

    xtitle = 'log (M / M_{\odot})'
    ytitle = '12 + log (O/H)'

    xrange = massrange
    yrange = ohrange
    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xthick=postthick1, ythick=postthick1, $
      charthick=postthick2, charsize=2.0, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle

    peg = peg1 & taucolor = taucolors[0] & tauline = taulines[0]
    x = alog10(peg.mstar*mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=0.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg3 & taucolor = taucolors[1] & tauline = taulines[1]
    x = alog10(peg.mstar*mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    peg = peg5 & taucolor = taucolors[2] & tauline = taulines[2]
    x = alog10(peg.mstar*mgalaxy) & y = peg.log12oh & age = peg.age+tform*1D3
    djs_oplot, x, y, color=fsc_color(taucolor,100), line=tauline, thick=postthick3
    xyouts, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), $
      string(labelz,format='(F3.1)'), charsize=2.0, align=1.0, charthick=postthick2
    plots, interpol(x,age,labelage*1D3), interpol(y,age,labelage*1D3), psym=8

    legend, textoidl('z_{f} = '+string(zform,format='(F3.1)')), /right, /bottom, $
      box=0, charsize=2.0, charthick=postthick2    
    legend, textoidl(['\tau'+tau_suffix+' = 1','\tau'+tau_suffix+' = 3','\tau'+tau_suffix+' = \infty']), $
      /left, /top, box=0, charsize=1.8, charthick=postthick2, spacing=2.2, $
      thick=postthick3, color=fsc_color(taucolors,[100,101,102]), line=taulines, $
      textcolor=fsc_color(taucolors,[100,101,102])

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

    if keyword_set(postscript) then dfpsclose

return
end
    
