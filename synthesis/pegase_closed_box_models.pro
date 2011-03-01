pro pegase_closed_box_models, postscript=postscript, encapsulated=encapsulated, paper=paper
; assumptions: (1) instantaneous recycling; (2) spatially and
; temporally invariant IMF; (3) constant yield

    if keyword_set(paper) then begin
       postscript = 1L
       encapsulated = 1L
    endif
    
;   paperpath = cwd()
    paperpath = ages_path(/papers)+'mz/FIG_MZ/'

    red, h100=0.7, omega0=0.3, omegal=0.7
    tuniverse = getage(0.0)
    
    alpha = 0.7          ; from stellar evolution
    mgalaxy = 1D0        ; final galaxy mass [M_sun]
    ginitial = 1D0       ; initial gas mass
    mass2number = 11.728 ; conversion
    yield = 0.0126       ; yield (Asplund et al. 2004)
    oyield = 7.4D-3      ; oxygen yield from Meynet & Maeder (2002)

; compute the full temporal evolution of the models on this time grid     
    
    maxtime = 13.5 & mintime = 1E-2 & dlogtime = 0.01
    logtime = findgen((alog10(maxtime)-alog10(mintime))/dlogtime+1)*dlogtime+alog10(mintime)
    time = 10.0^logtime
;   maxtime = 13.5 & dtime = 0.01 & mintime = 1E-3
;   time = findgen((maxtime-mintime)/dtime+1)*dtime+mintime ; [Gyr]
    redshift = getredshift(time)

; consider two formation redshifts, z1=7 and z2=2

    zform1 = 7.0 & zform2 = 2.0 
    tform1 = getage(zform1) & age1 = tuniverse-tform1
    tform2 = getage(zform2) & age2 = tuniverse-tform2

    redshift1 = getredshift(time+tform1)
    redshift2 = getredshift(time+tform2)
    
    splog, 'zform1 = '+string(zform1,format='(F3.1)')+', age1 = '+string(age1,format='(F5.2)')+' Gyr, '+$
      'tform1 = '+string(tform1,format='(F4.2)')+' Gyr'
    splog, 'zform2 = '+string(zform2,format='(F3.1)')+', age2 = '+string(age2,format='(F5.2)')+' Gyr, '+$
      'tform2 = '+string(tform2,format='(F4.2)')+' Gyr'

    tau = [2.0,5.0,10.0,100.0] ; [Gyr]
;   tau = [2.0,3.0,5.0,10.0,100.0] ; [Gyr]
;   tau = [0.5,1.0,2.0,3.0,5.0,8.0,10.0,20.0,40.0,1E3] ; [Gyr]
    ntau = n_elements(tau)
;   color = ['','blue','dark green','red','dark cyan']
    color = ['','dark green','red','navy']
    line = [0,2,3,4]
    if (ntau ne n_elements(color)) or (ntau ne n_elements(line)) then begin
       splog, 'Dimensions must match!'
       return
    endif

    gasrange = [-1.5,0.3]
    fracrange = [0.0,1.0]
    Zrange = [0.0,0.015]
    ohrange = [7.0,9.1]

;   yrange = gasrange
;   yrange = fracrange
;   yrange = Zrange
    yrange = ohrange
    xrange = [0.0,1.5]
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.1,0.9], ymargin=[0.9,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    if keyword_set(postscript) then begin
       postthick = 5.0
       postthick2 = 8.0
       postthick3 = 4.0
    endif else begin
       postthick = 2.0
       postthick2 = 2.0
       postthick3 = 1.5
    endelse
    
; ##################################################
; closed box (A=E=0); exponentially declining SFR
; ##################################################

    psname = 'closed_box_sfr_tau'
    if keyword_set(encapsulated) then suffix = '.eps' else suffix = '.ps'
    if keyword_set(postscript) then dfpsplot, paperpath+psname+suffix, /color, $
      xsize=8.5, ysize=8.5, encapsulated=encapsulated

    for i = 0L, ntau-1L do begin

; formation redshift = infinity
       
       gas = alog(1.0 - (alpha*mgalaxy/ginitial)*(1.0-exp(-time/tau[i])))
       frac = 1.0 - (alpha*mgalaxy/ginitial)*(1.0-exp(-time/tau[i]))
       Z = yield*(time/tau[i] - alog((alpha*mgalaxy/ginitial)*(1.0-exp(time/tau[i]))+exp(time/tau[i])))
       oh = 12.0+alog10(Z*oyield/yield/mass2number)

       oh_predict = 12.0 + alog10((oyield/mass2number)*alog(1.0/frac))
;      plot, oh_predict, oh-oh_predict, ysty=3 & cc = get_kbrd(1)
       
       if (i eq 0L) then begin
          djs_plot, [0], [0], /nodata, xsty=8, ysty=11, xrange=xrange, $
            yrange=yrange, xthick=postthick, ythick=postthick, charsize=1.8, $
            charthick=postthick, xtitle='Redshift', ytitle='12 + log (O/H)', position=pos
          xrange2 = getage(0.0)-getage(xrange)
          axis, /xaxis, xthick=postthick, xsty=1, xrange=xrange2, xtitle='Lookback time [Gyr]', $
            charsize=1.8, charthick=postthick
          yrange2 = exp(-(mass2number/oyield)*10.0^(!y.crange-12.0))
          ytitle2 = 'Gas Fraction' ; textoidl('ln (1/\mu)')
;         srt = sort(oh) & linterp, oh[srt], alog(1.0/frac[srt]), !y.crange, yrange2
          axis, /yaxis, ythick=postthick, ysty=3, yrange=yrange2, $
            ytitle=ytitle2, charsize=1.8, charthick=postthick
       endif

;      djs_oplot, redshift, oh, line=i, thick=postthick2, color=color[i]
       djs_oplot, redshift1, oh, line=line[i], thick=postthick2, color=color[i]
       djs_oplot, redshift2, oh, line=line[i], thick=postthick2, color=color[i]
       print, tau[i], interpol(oh,redshift1,0.0)-interpol(oh,redshift1,1.0), $
         interpol(oh,redshift2,0.0)-interpol(oh,redshift2,1.0);, minmax(frac), minmax(oh)

       zpos = 1.1
       case i of
          0L: xyouts, zpos, 8.68, textoidl('\tau = '+string(tau[i],format='(I0)')), $
            align=0.5, charsize=1.5, charthick=postthick, /data
          1L: xyouts, zpos, 8.40, textoidl('\tau = '+string(tau[i],format='(I0)')), $
            align=0.5, charsize=1.5, charthick=postthick, /data
          2L: xyouts, zpos, 8.10, textoidl('\tau = '+string(tau[i],format='(I0)')), $
            align=0.5, charsize=1.5, charthick=postthick, /data
          3L: xyouts, zpos, 7.15, textoidl('\tau = '+string(tau[i],format='(I0)')), $
            align=0.5, charsize=1.5, charthick=postthick, /data
          else: 
       endcase

;      xyouts, 0.1, 8.4, textoidl('\psi(t) = \psi_{0} e^{-t/\tau}'), align=0.0, $
;        charsize=2.0, charthick=postthick3, /data
       xyouts, 0.2, 8.0, textoidl('\psi(t) \propto e^{-t/\tau}'), align=0.0, $
         charsize=2.0, charthick=postthick3, /data
       xyouts, 0.2, 7.8, 'Closed box', align=0.0, charsize=2.0, charthick=postthick3, /data
       
;      case i of
;         0L: xyouts, 0.9, 9.35, textoidl('\tau = '+string(tau[i],format='(I0)')), $
;           align=0.5, charsize=1.5, charthick=postthick, /data
;         1L: xyouts, 0.6, 9.15, textoidl('\tau = '+string(tau[i],format='(I0)')), $
;           align=0.5, charsize=1.5, charthick=postthick, /data
;         2L: xyouts, 0.4, 9.00, textoidl('\tau = '+string(tau[i],format='(I0)')), $
;           align=0.5, charsize=1.5, charthick=postthick, /data
;         3L: xyouts, 0.2, 8.15, textoidl('\tau = '+string(tau[i],format='(I0)')), $
;           align=0.5, charsize=1.5, charthick=postthick, /data
;         else: 
;      endcase

;      case i of
;         0L: xyouts, 0.9, 9.35, textoidl('\tau = '+string(tau[i],format='(I0)')), $
;           align=0.5, charsize=1.5, charthick=postthick, /data
;         1L: xyouts, 0.6, 9.15, string(tau[i],format='(I0)'), align=0.5, charsize=1.5, charthick=postthick, /data
;         2L: xyouts, 0.4, 9.00, string(tau[i],format='(I0)'), align=0.5, charsize=1.5, charthick=postthick, /data
;         3L: xyouts, 0.2, 8.15, string(tau[i],format='(I0)'), align=0.5, charsize=1.5, charthick=postthick, /data
;         else: 
;      endcase

;       if (i eq 0L) then begin
;          srt = sort(oh) & linterp, oh[srt], alog(1.0/frac[srt]), !y.crange, yrange2
;;         yrange2 = interpol(alog(1.0/frac),oh,!y.crange)
;          axis, /yaxis, ythick=postthick, ysty=3, yrange=yrange2, $
;            xtitle='ln (1/\mu)', charsize=1.8, charthick=postthick
;       endif
;          
;       djs_plot, redshift, alog(1.0/frac), line=2, thick=2.0, /noerase, $
;         xsty=4, ysty=4, xrange=xrange, yrange=yrange2, position=pos

;      djs_oplot, time, Z, line=0, thick=2.0
;      djs_oplot, time, gas, line=0, thick=2.0
;      print, tau[i] & cc = get_kbrd(1)

; verify the closed box model:       

;      plot, alog(1.0/frac), z, thick=2 & print, linfit(alog(1.0/frac),z)
       
    endfor

    if keyword_set(postscript) then dfpsclose

stop
    

return
end
