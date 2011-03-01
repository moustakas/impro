pro render_gandalf_lineplot, wave, flux_noc, linefit, specdata, $
  noerase=noerase, xrange=xrange, scale=scale, pos=pos, $
  linelabel=linelabel

    get_element, wave, xrange, xx
    max_norej = max(linefit[xx[0]:xx[1]])
;   max_norej = max(flux_noc[xx[0]:xx[1]])
    max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
    sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)

    yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
      (max_norej*1.25)>(3.0*sigma_rej)]
;   yrange = (im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)>$
;     max(linefit[xx[0]:xx[1]]))*[-0.12,1.2]
    
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos, ytitle='', noerase=noerase, $
      xtickinterval=25
    djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
    djs_oplot, wave, scale*linefit, ps=10, color='red', thick=2.5
    oplot_gandalf_linewave, specdata
    im_legend, linelabel, /left, /top, box=0, charsize=1.2

return
end

pro oplot_gandalf_linewave, specdata
; overplot the wavelength of the line, adjusting for the possible
; differences in the emission/absorption-line redshifts
    for jj = 0, n_elements(specdata.linename)-1 do begin
       line = strupcase(strtrim(specdata.linename[jj],2))
       linewave = specdata.(tag_indx(specdata,line+'_WAVE'))
       linez = (specdata.(tag_indx(specdata,line+'_LINEZ')))[0]
       if (linez eq 0.0) then factor = 1.0 else $
         factor = (1.0+linez)/(1.0+specdata.zabs)
       djs_oplot, factor*linewave*[1,1], !y.crange, line=5, color='blue'
    endfor
    djs_oplot, !x.crange, [0,0], line=0, color='red'
return
end

;;pro render_gandalf_lineplot, wave, flux_noc, linefit, specdata, $
;;  noerase=noerase, xrange=xrange, scale=scale, pos=pos, $
;;  linelabel=linelabel
;;
;;    get_element, wave, xrange, xx
;;    max_norej = max(linefit[xx[0]:xx[1]])
;;;   max_norej = max(flux_noc[xx[0]:xx[1]])
;;    max_rej = im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)
;;    sigma_rej = djsig(flux_noc[xx[0]:xx[1]],sigrej=3.0)
;;
;;    yrange = [(max_rej*(-0.12))<(-2.5*sigma_rej),$
;;      (max_norej*1.25)>(3.0*sigma_rej)]
;;;   yrange = (im_max(flux_noc[xx[0]:xx[1]],sigrej=3.0)>$
;;;     max(linefit[xx[0]:xx[1]]))*[-0.12,1.2]
;;    
;;    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
;;      yrange=scale*yrange, position=pos, ytitle='', noerase=noerase, $
;;      xtickinterval=25
;;    djs_oplot, wave, scale*flux_noc, ps=10, color='grey', thick=3
;;    djs_oplot, wave, scale*linefit, ps=10, color='red', thick=2.5
;;    oplot_linewave, specdata
;;    im_legend, linelabel, /left, /top, box=0, charsize=1.2
;;
;;return
;;end
;;
;;pro oplot_linewave, specdata
;;    for jj = 0, n_elements(specdata.linename)-1 do djs_oplot, $
;;      specdata.(tag_indx(specdata,strupcase(strtrim(specdata.linename[jj],2))+$
;;      '_WAVE'))*[1,1], !y.crange, line=5, color='blue'
;;    djs_oplot, !x.crange, [0,0], line=0, color='red'
;;return
;;end

pro qaplot_gandalf, specdata, specfit, pos1, pos2
; make the actual QAplot - see below for the main routine
    
    light = 2.99792458D5 ; speed of light [km/s]

; define some convenient internal variables
    galaxy = 'ages '+string(specdata.pass,format='(I3.3)')+'/'+$
      string(specdata.aper,format='(I3.3)')
    
    good = where(specfit.wave gt 0.0)
    wave = exp(specfit.wave[good])
    flux = specfit.flux[good]
    linefit = specfit.linefit[good]
    continuum = specfit.continuum[good]
    smooth_continuum = specfit.smooth_continuum[good]
    flux_noc = flux - continuum - smooth_continuum

    xtitle1 = 'Rest Wavelength (\AA)'
    power = ceil(abs(alog10(median(flux))))
    scale = 10.0^power
    ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    
; ---------------------------------------------------------------------------    
; page1 - continuum plot

; first third    
    xrange = [min(wave),(max(wave)-min(wave))/3.0+min(wave)]
;   xrange = [min(wave),mean(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=xrange, $
      yrange=scale*yrange, position=pos1[*,0], ytitle='', $
      xtitle='', title=galaxy
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5

    im_legend, /right, /top, box=0, margin=0, charsize=1.2, $
      ['\chi^2_{\nu}='+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2),$
      'E(B-V)='+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2),$;+'\pm'+$
;     strtrim(string(specdata.continuum_ebv_err,format='(F12.3)'),2),$
      '\sigma='+strtrim(string(specdata.vdisp,format='(F12.1)'),2)+' km s^{-1}',$
;     '\sigma='+strtrim(string(specdata.vdisp,format='(F12.1)'),2)+'\pm'+$
;     strtrim(string(specdata.vdisp_err,format='(F12.1)'),2)+' km s^{-1}']
      'S/N='+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)]

    im_legend, /left, /top, box=0, margin=0, charsize=1.2, $
      ['z_{ages}='+strtrim(string(specdata.z,format='(F12.5)'),2),$
      'z_{abs}='+strtrim(string(specdata.zabs,format='(F12.5)'),2)]
;     'z_{Balmer}='+strtrim(string(specdata.zline_balmer,format='(F12.5)'),2),$
;     'z_{forbid}='+strtrim(string(specdata.zline_forbidden,format='(F12.5)'),2)]

; second third
    xrange = [(max(wave)-min(wave))/3.0+min(wave),(max(wave)-min(wave))*2.0/3.0+min(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, $
      xrange=xrange, yrange=scale*yrange, position=pos1[*,1], $
      ytitle=ytitle1, xtitle=''
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5
    
; third third    
    xrange = [(max(wave)-min(wave))*2.0/3.0+min(wave),max(wave)]
;   xrange = [mean(wave),max(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, $
      xrange=xrange, yrange=scale*yrange, position=pos1[*,2], $
      ytitle='', xtitle=xtitle1
    djs_oplot, wave, scale*flux, ps=10, color='grey'
    djs_oplot, wave, scale*continuum, ps=10, color='red', thick=2.5
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      color='blue', thick=2.5
    djs_oplot, wave, scale*smooth_continuum, ps=10, color='dark green', thick=2.5
    
; ---------------------------------------------------------------------------    
; page2 - emission-line plot

    npanel = 8
    linelabel = ['[Ne V]','[O II]','H\delta','H\gamma',$
      'H\beta','[O III]','H\alpha+[N II]','[S II]']
    xrange = fltarr(2,npanel)
    xrange[*,0] = 3426.0
    xrange[*,1] = 3727.0
    xrange[*,2] = 4101.0
    xrange[*,3] = 4340.0
    xrange[*,4] = 4861.0
    xrange[*,5] = [4959.0,5007.0]
    xrange[*,6] = [6548.0,6584.0]
    xrange[*,7] = [6716.0,6731.0]
    xrange = xrange+rebin(30.0*[-1,1],2,npanel)

    for ii = 0, npanel-1 do begin
       render_gandalf_lineplot, wave, flux_noc, linefit, $
         specdata, xrange=xrange[*,ii], scale=scale, $
         linelabel=linelabel[ii], pos=pos2[*,ii], noerase=(ii gt 0)
    endfor
    
return
end
    
pro qaplot_gandalf_specfit, pass1, test=test, doplot=doplot
; jm09nov16ucsd - build a QAplot from the AGES_GANDALF_SPECFIT output

; path names and emission-line file name
    version = ages_version(/ppxf_specfit)
    base_specfitpath = ages_path(/ppxf)
    specfitpath = base_specfitpath+'fluxed/tweak/'+version+'/' 

    if keyword_set(test) then suffix = '_test' else suffix = ''
    if (n_elements(pass1) eq 0) then begin
       specdatafiles = file_search(specfitpath+'specdata_???'+$
         suffix+'.fits.gz',count=nall)
       specfitfiles = repstr(specdatafiles,'specdata','specfit')
    endif else begin
       pass = string(pass1,format='(I3.3)')
       specdatafiles = file_search(specfitpath+'specdata_'+pass+$
         suffix+'.fits.gz',count=nall)
       specfitfiles = repstr(specdatafiles,'specdata','specfit')
    endelse
    psfiles = specfitpath+'qaplot_'+file_basename(repstr(repstr($
      specdatafiles,'.fits.gz','.ps'),'specdata_',''))

    im_plotconfig, 19, pos2, xspace=0.5, yspace=[0.5,0.5,0.5], $
       width=[3.25,3.25], height=[2.0,2.0,2.0,2.0]
    for ii = 0, nall-1 do begin
       im_plotconfig, 4, pos1, psfile=psfiles[ii], $
         yspace=[0.5,0.5], charsize=1.2, ymargin=[0.6,1.0], $
         thick=2.0, height=replicate(2.8,3)
       
       specdata = mrdfits(specdatafiles[ii],1)
       specfit = mrdfits(specfitfiles[ii],1)

       index = lindgen(n_elements(specdata))
;      index = where((specdata.nev_3426[0]/specdata.nev_3426[1] gt 3.0))
;      index = where((specdata.h_beta[0]/specdata.h_beta[1] gt 3.0))
;      index = where(3E5*abs(specdata.zabs-specdata.z) gt 200.0)

       if (index[0] ne -1) then begin
;         for jj = 0, 5 do begin
          for jj = 0, n_elements(index)-1 do begin
             print, format='("Object ",I0,"/",I0,A10,$)', jj+1, $
               n_elements(index), string(13b)
             qaplot_gandalf, specdata[index[jj]], $
               specfit[index[jj]], pos1, pos2
          endfor
       endif
       
       im_plotconfig, psfile=psfiles[ii], /gzip, /psclose
       spawn, 'rsync -auv '+psfiles[ii]+'.gz ~', /sh
    endfor

return
end    
