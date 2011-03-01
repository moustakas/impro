;+
; NAME:
;   QAPLOT_ISEDFIT_SFHGRID
;
; PURPOSE:
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;
; Copyright (C) 2010, John Moustakas
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

function init_info, nmodel, nage=nage
    info = {$
      ub: fltarr(nage),$
      bv: fltarr(nage)}
    info = replicate(info,nmodel)
return, info
end

pro qaplot_isedfit_sfhgrid, synthmodels=synthmodels, imf=imf, $
  sfhgrid=sfhgrid, redcurve=redcurve

    common com_info, info
    
    if (n_elements(synthmodels) eq 0L) then synthmodels = 'bc03'
    if (n_elements(imf) eq 0L) then imf = 'chab' ; 'salp'
    if (n_elements(sfhgrid) eq 0L) then sfhgrid = 2
    if (n_elements(redcurve) eq 0L) then redcurve = 0
    sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')

; reddening curves    
    case redcurve of
       0: redcurvestring = 'calzetti'
       1: redcurvestring = 'charlot'
       2: redcurvestring = 'odonnell'
       3: redcurvestring = 'smc'
       else: begin
          splog, 'Unrecognized reddening curve'
          return
       end
    endcase
    
; build the absolute path names    
    isedfit_sfhgridpath = getenv('ISEDFIT_SFHGRID_DIR')+'/'
    case synthmodels of
       'bc03': begin
          basemodelspath = isedfit_sfhgridpath+'basemodels/bc03/'
          sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')
          sfhgridpath = isedfit_sfhgridpath+sfhgridstring+$
            '/bc03/'+redcurvestring+'/'
          baseinfo = mrdfits(isedfit_sfhgridpath+'basemodels/info_bc03_'+$
            imf+'.fits.gz',1,/silent)
       end
       else: message, 'Stellar pops models not supported!'
    endcase

; build the info structure we need
    infofile = sfhgridpath+imf+'_qainfo.fits'

    chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits.gz'
    if (file_test(chunkinfofile,/regular) eq 0L) then $
      message, 'No SFH grid information file found!'
    chunkinfo = mrdfits(chunkinfofile,1,/silent)
    chunkfiles = sfhgridpath+strtrim(chunkinfo.chunkfiles,2)
    nchunk = n_elements(chunkfiles)

    if (n_elements(info) eq 0L) then begin
       for ichunk = 0L, nchunk-1L do begin
          splog, 'Reading '+chunkfiles[ichunk]
          info1 = mrdfits(chunkfiles[ichunk],1)
          info1 = struct_trimtags(info1,except=['wave','flux'])
          if (n_elements(info) eq 0L) then info = info1 else $
            info = [temporary(info),info1]
       endfor         
    endif
    nmodel = n_elements(info)
    
; color-color diagram of the full range of models
;   nuvmr = reform(info.abmag[2,*]-info.abmag[3,*])
;   rmi = reform(info.abmag[3,*]-info.abmag[4,*])
    nuvmr = reform(info.abmag[1,*]-info.abmag[5,*])
    rmi = reform(info.abmag[5,*]-info.abmag[6,*])
    hasburst = where(info.nburst gt 0,nhasburst,$
      comp=issmooth,ncomp=nissmooth)
       
    psfile = sfhgridpath+imf+'_qainfo.ps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=6.8
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='R - I', ytitle='NUV - R', xrange=[-0.5,0.8], $
      yrange=[-2.2,8.0]
    if (nhasburst ne 0L) then djs_oplot, rmi[*,hasburst], $
      nuvmr[*,hasburst], psym=symcat(9), symsize=0.2, color='red'
    if (nissmooth ne 0L) then djs_oplot, rmi[*,issmooth], $
      nuvmr[*,issmooth], psym=symcat(16), symsize=0.2, color='blue'
    im_plotconfig, psfile=psfile, /psclose, /gzip

    
stop    
    
; color-color and color-magnitude diagrams at various redshift bins
; assuming a fixed formation redshift
;   zf = 2.0
;   zbins = [0.1,0.3,0.5,0.7,0.9]
    zf = 3.0
    zbins = [0.4,0.9]
    nz = n_elements(zbins)
    psfile = sfhgridpath+imf+'_qainfo.ps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=6.8
    for iz = 0, nz-1 do begin
       thisage = getage(zbins[iz])-getage(zf)
       splog, zbins[iz], thisage
       aindx = findex(info[0].age,thisage)
       hasburst = where(info.nburst gt 0,nhasburst,$
         comp=issmooth,ncomp=nissmooth)

       nuvmr = interpolate(transpose(reform(info.abmag[1,*] - info.abmag[5,*])),aindx)
       rmi = interpolate(transpose(reform(info.abmag[5,*] - info.abmag[6,*])),aindx)
       
       djs_plot, [0], [0], /nodata, position=pos, $
         xsty=1, ysty=1, xtitle='R - I', ytitle='NUV - R', $
         xrange=[-0.5,0.8], yrange=[-2.2,8.0], title='z='+$
         string(zbins[iz],format='(F4.2)')
       if (nhasburst ne 0L) then djs_oplot, rmi[hasburst], $
         nuvmr[hasburst], psym=symcat(9,thick=5), color='red'
       if (nissmooth ne 0L) then djs_oplot, rmi[issmooth], $
         nuvmr[issmooth], psym=symcat(16,thick=5), color='blue'
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
stop    
    
; read all the models    
    nmodel = n_elements(chunk)
    nage = n_elements(chunk[0].modelage)
    bigage = chunk.modelage
    
    bigtau = rebin(reform(chunk.tau,1,nmodel),nage,nmodel)
    bigebv = rebin(reform(chunk.ebv,1,nmodel),nage,nmodel)
    bigsfr = bigage*0.0
    notssp = where(bigtau ne 0.0)
    bigsfr[notssp] = exp(-bigage[notssp]/bigtau[notssp])
    
    plot, bigebv, (bigtau>1E-3)/bigage, psym=6, xsty=3, ysty=3, /ylog

    ebvgrid = im_array(0.0,1.0,0.02)
    tau = range(1E-1,100,100,/log)
    age = range(1E-2,12.0,100,/log)
    ratio = age/tau

    ratio = bigtau/bigage
;   ratio = range(1E-3,1E4,500,/log)
    plot, ratio, exp(-ratio/0.25), xsty=3, ysty=3, $
      yrange=[0,1], ps=6;, /xlog, xrange=[0.01,100]
    djs_oplot, 0.25*[1,1], !y.crange
    
    
    
stop    
    
    ugriz = fltarr(5,nage,nmodel)
    for ii = 0L, nmodel-L do begin
       k_projection_table, rmatrix, chunk.flux, k_lambda_to_edges(chunk[0].wave), $
         0.0, sdss_filterlist()
       
       ri_tau = -2.5*alog10(rmatrix[0,*,1]/rmatrix[0,*,2])
       gr_tau = -2.5*alog10(rmatrix[0,*,0]/rmatrix[0,*,1])

       
    endfor

; diagnostic plots                      
                if keyword_set(plotcolors) then begin
                   oo = outinfo[indx1[jj]]
                   k_projection_table, rmatrix, tauflux, k_lambda_to_edges(taufits.wave), $
                     0.0, 'sdss_'+['g0','r0','i0']+'.par'
                   ri_tau = -2.5*alog10(rmatrix[0,*,1]/rmatrix[0,*,2])
                   gr_tau = -2.5*alog10(rmatrix[0,*,0]/rmatrix[0,*,1])

                   k_projection_table, rmatrix, outflux, k_lambda_to_edges(taufits.wave), $
                     0.0, 'sdss_'+['g0','r0','i0']+'.par'
                   ri = -2.5*alog10(rmatrix[0,*,1]/rmatrix[0,*,2])
                   gr = -2.5*alog10(rmatrix[0,*,0]/rmatrix[0,*,1])

                   djs_plot, ri_tau, gr_tau, ps=-4, xtitle='r-i', ytitle='g-r', color='blue'
                   djs_oplot, ri, gr, ps=-6, color='green'
                   get_element, outage, tburst0, kk
                   plots, ri[kk], gr[kk], psym=7, color=djs_icolor('red'), sym=2
                   legend, ['TAU = '+strtrim(oo.tau,2),'TFORM = '+strtrim(oo.tform,2),$
                     'TBURST = '+strtrim(oo.tburst,2),'FBURST = '+strtrim(oo.fburst,2),$
                     'DTAUBURST = '+strtrim(oo.dtauburst,2)], $
                     /left, /top, box=0, charsize=1.6
                   cc = get_kbrd(1)
                endif

                if keyword_set(plotsfr) then begin
                   oo = outinfo[indx1[jj]]
                   if (oo.tau gt 0.0) then begin
                      if (ichunk eq 0L) then window, 2 else wset, 2
                      djs_plot, iage, tausfr, ps=-4, xsty=3, color='red', $
                        yr=minmax(outsfr), xrange=[0.0,maxage], /ylog, $
                        xtitle='Time since Formation', ytitle='SFR'
                      djs_oplot, iage, burstsfr, color='blue'
                      djs_oplot, outage, outsfr, ps=-4
                      oo = outinfo[indx1[jj]]
                      legend, ['TAU = '+strtrim(oo.tau,2),'TFORM = '+strtrim(oo.tform,2),$
                        'TBURST = '+strtrim(oo.tburst,2),'FBURST = '+strtrim(oo.fburst,2),$
                        'DTAUBURST = '+strtrim(oo.dtauburst,2)], $
                        /right, /top, box=0, charsize=1.6
                      if keyword_set(plotspec) then begin
                         if (ichunk eq 0L) then window, 0 else wset, 0
                         for kk = 0L, n_elements(outage)-1L do begin
                            djs_plot, taufits.wave, outflux[*,kk], xrange=[3000,9000]
                            legend, ['AGE = '+strtrim(outage[kk],2),$
                              'AGE_UNIVERSE = '+strtrim(oo.tform+outage[kk],2)], $
                              /right, /top, box=0, charsize=1.6
                            cc = get_kbrd(1)
                         endfor
                      endif else cc = get_kbrd(1)
                   endif
                endif


; diagnostic plots                      
                if keyword_set(plotcolors) then begin
                   oo = outinfo[indx1[jj]]
                   k_projection_table, rmatrix, outflux, k_lambda_to_edges(taufits.wave), $
                     0.0, 'sdss_'+['g0','r0','i0']+'.par'
                   ri = -2.5*alog10(rmatrix[0,*,1]/rmatrix[0,*,2])
                   gr = -2.5*alog10(rmatrix[0,*,0]/rmatrix[0,*,1])

                   djs_plot, ri, gr, ps=-4, xtitle='r-i', ytitle='g-r', color='blue'
                   legend, ['TAU = '+strtrim(oo.tau,2),'TFORM = '+strtrim(oo.tform,2)], $
                     /left, /top, box=0, charsize=1.6
                   cc = get_kbrd(1)
                endif

    
    
    
stop    

return
end
;;    
;;       filterlist = 'bessell_'+['U','B','V']+'.par'
;;       nfilt = n_elements(filterlist)
;;       for ichunk = 0L, nchunk-1L do begin
;;          splog, 'Reading '+chunkfiles[ichunk]
;;          chunk1 = mrdfits(chunkfiles[ichunk],1,range=[0,5])
;;          nmodel = n_elements(chunk1)
;;          npix = n_elements(chunk1[0].wave)
;;          info = init_info(nmodel,nage=params[0].nage)
;;
;;          if (ichunk eq 0L) then edgewave = k_lambda_to_edges(chunk1[0].wave)
;;          k_projection_table, rmatrix, reform(chunk1.flux,npix,nmodel*params[0].nage), $
;;            /silent, edgewave, zvals, filterlist, zmin=0.0, zmax=0.0, nz=1
;;          maggies = transpose(reform(rmatrix,nmodel,params[0].nage,nfilt),[1,0,2])
;;          info.ub = -2.5*alog10(maggies[*,*,0]/maggies[*,*,1])
;;          info.bv = -2.5*alog10(maggies[*,*,1]/maggies[*,*,2])
;;
;;          
;;stop
;;; technically we don't have to loop, but it makes the problem
;;; just a bit more manageable
;;;         for iage = 0, params[0].nage-1 do begin
;;;            k_projection_table, rmatrix, chunk1.flux[*,iage], /silent, $
;;;              edgewave, zvals, filterlist, zmin=0.0, zmax=0.0, nz=1
;;;         endfor
;;       endfor
;;       
