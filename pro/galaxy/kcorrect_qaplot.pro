;+
; NAME:
;   KCORRECT_QAPLOT
;
; PURPOSE:
;   Build a QAplot of the K-correct SED fitting.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is generally called from within IM_KCORRECT, but it
;   can also be called indepedently.  See IM_KCORRECT for an example
;   of the input syntax.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jul 06, NYU
;   jm10jan24ucsd - documented and generalized
;
; Copyright (C) 2009-2010, John Moustakas
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

function kcorrect_restore, info, vname=vname
; internal function to rebuild the best-fitting SED in the observed 
; frame  

    if (n_elements(vname) eq 0) then vname = 'default.nolines'
    ngal = n_elements(info)
    
    light = 2.99792458D18       ; speed of light [A/s]

    k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
      lfile=lfile, vpath=vpath, vname=vname
    wave = k_lambda_to_centers(lambda)
    model = {wave: wave, flux: wave*0.0}
    model = replicate(temporary(model),ngal)
    model.flux = vmatrix#info.k_coeffs
    for igal = 0L, ngal-1L do begin
       model[igal].wave = wave*(1.0+info[igal].k_zobj) ; [A]
       model[igal].flux = model[igal].flux/(1.0+info[igal].k_zobj)
    endfor

    model.flux = model.flux*model.wave^2/light ; observed frame
    model.flux = -2.5*alog10(model.flux>1D-50)-48.6 ; [AB mag]
    
return, model
end

pro kcorrect_qaplot, info, psfile=psfile, vname=vname, $
  in_filterlist=in_filterlist, clobber=clobber, noplot=noplot
    
    ngal = n_elements(info)
    if (ngal eq 0L) then begin
       doc_library, 'kcorrect_qaplot'
       return
    endif

; restore the filters and best-fitting models
    if (n_elements(in_filterlist) eq 0) then $
      in_filterlist = strtrim(info[0].filterlist,2)
    in_filtinfo = im_filterspecs(filterlist=in_filterlist)
    in_nfilt = n_elements(in_filterlist)

    if tag_exist(info,'out_filterlist') then begin
       out_filterlist = strtrim(info[0].out_filterlist,2)
       out_filtinfo = im_filterspecs(filterlist=out_filterlist)
       out_nfilt = n_elements(out_filterlist)
    endif

    model = kcorrect_restore(info,vname=vname)
    if keyword_set(noplot) then return

    if tag_exist(info,'galaxy') then galaxy = info.galaxy else begin
       fmt = '(I'+string(5L,format='(I0)')+'.'+string(5L,format='(I0)')+')'
       galaxy = 'Galaxy_'+string(lindgen(ngal),format=fmt)
    endelse

; generate the QA-plot
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    if (n_elements(psfile) ne 0) then begin
       if file_test(psfile+'.gz',/regular) and $
         (keyword_set(clobber) eq 0) then begin
          splog, 'Output file '+psfile+' exists; use /CLOBBER'
          return
       endif
    endif

    im_plotconfig, 8, pos, psfile=psfile, ymargin=[1.1,1.1]
    for igal = 0L, ngal-1L do begin
       if ((igal mod 10) eq 0) then print, igal+1L, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'
       
; redshift and best-fit model
       zobj = info[igal].k_zobj
       wave = model[igal].wave    ; [A]
       flux = model[igal].flux    ; [AB]

; make the plot       
       in_xrange1 = [min(in_filtinfo.weff-1.3*in_filtinfo.fwhm),$
         max(in_filtinfo.weff+2.1*in_filtinfo.fwhm)]
       if tag_exist(info,'out_filterlist') then begin
          out_xrange1 = [min(out_filtinfo.weff-1.3*out_filtinfo.fwhm),$
            max(out_filtinfo.weff+2.1*out_filtinfo.fwhm)]
          xrange1 = [in_xrange1[0]<out_xrange1[0],in_xrange1[1]>out_xrange1[1]]
       endif else xrange1 = in_xrange1
       
       xrange1[0] = xrange1[0]>(500.0*(1+zobj))
       xrange1[1] = xrange1[1]<max(wave)

       if (info[igal].k_chi2 eq -999.0) then begin
          xrange2 = xrange1/(1.0+zobj)
          djs_plot, [0], [0], /nodata, xsty=9, ysty=1, xrange=xrange1, $
            yrange=yrange, /xlog, xtitle=xtitle1, ytitle=ytitle1, $
            ytickname=replicate(' ',10), position=pos        ;, xtickinterval=1000
          axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog
          legend, ['No mass estimate available'], /left, /top, $
            box=0, spacing=1.5, charsize=1.6
          
          label = [strtrim(galaxy[igal],2),'z = '+string(zobj,format='(F6.4)')]
          legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.6
       endif else begin
          get_element, wave, xrange1, xx

          bestmab = -2.5*alog10(info[igal].k_bestmaggies)

          yrange = fltarr(2)
          yrange[0] = (max(bestmab)>max(flux[xx[0]:xx[1]]))*1.05
          yrange[1] = (min(bestmab)<min(flux[xx[0]:xx[1]]))*0.95

          djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
            xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
            position=pos        ;, xtickinterval=1000
          axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog

          djs_oplot, wave, flux, line=0, color='grey'
          djs_oplot, in_filtinfo.weff, bestmab, $
            psym=symcat(6,thick=6), symsize=2.5

; overplot the data; distinguish between three different cases, based
; on the input photometry
          used = where((info[igal].k_maggies gt 0.0) and $ ; used in the fitting
            (info[igal].k_ivarmaggies gt 0.0),nused)
          notused = where((info[igal].k_maggies gt 0.0) and $ ; not used in the fitting
            (info[igal].k_ivarmaggies eq 0.0),nnotused)
          nodata = where((info[igal].k_maggies eq 0.0) and $ ; no measurement
            (info[igal].k_ivarmaggies eq 0.0),nnodata)

          if (nused ne 0L) then begin
             mab = maggies2mag(info[igal].k_maggies[used],$
               ivar=info[igal].k_ivarmaggies[used],magerr=mab_err)
             oploterror, in_filtinfo[used].weff, mab, in_filtinfo[used].fwhm/2.0, $
               mab_err, psym=symcat(16), symsize=2.0, color=djs_icolor('dark green'), $
               errcolor=djs_icolor('dark green'), errthick=!p.thick
          endif

          if (nnotused ne 0L) then begin
             mab = maggies2mag(info[igal].k_maggies[notused])
             oploterror, in_filtinfo[notused].weff, mab, in_filtinfo[notused].fwhm/2.0, $
               mab*0.0, psym=symcat(4,thick=6.0), symsize=3.0, color=djs_icolor('red'), $
               errcolor=djs_icolor('red'), errthick=!p.thick
          endif

; overlay the observed *output* photometry
          if tag_exist(info,'outmaggies_obs') then begin
             outmab = maggies2mag(info[igal].outmaggies_obs,$
               ivar=info[igal].outivarmaggies_obs,magerr=outmab_err)
             oploterror, out_filtinfo.weff, outmab, out_filtinfo.fwhm/2.0, $
               outmab_err, psym=symcat(14), symsize=2.5, color=djs_icolor('orange'), $
               errcolor=djs_icolor('orange'), errthick=!p.thick
          endif
          
; overlay the synthesized output photometry
          if tag_exist(info,'synth_outmaggies_obs') then begin
             synth_outmab = -2.5*alog10(info[igal].synth_outmaggies_obs)
             djs_oplot, out_filtinfo.weff, synth_outmab, $
               psym=symcat(5,thick=6), symsize=2.5, color='blue'
          endif
          
; legend
;         label = ['log (M/M'+sunsymbol()+') = '+$
;           strtrim(string(info[igal].mass,format='(F12.2)'),2)]
;         legend, textoidl(label), /left, /top, box=0, spacing=1.7, charsize=1.6
          label = textoidl([strtrim(repstr(galaxy[igal],'_',' '),2),$
            'z = '+string(zobj,format='(F6.4)'),'\chi^{2} = '+$
            strtrim(string(abs(info[igal].k_chi2),format='(F12.2)'),2),$
            'log (M/M'+sunsymbol()+') = '+$
            strtrim(string(info[igal].k_mass,format='(F12.2)'),2)])
          legend, label, /left, /top, box=0, spacing=1.5, charsize=1.4

;         mylabel = ['SDSS/ugriz','Model/ugriz','NDWFS/BwRI','Model/BwRI']
;         psym = [16,6,14,5]
;         color = ['dark green','','orange','blue']
;         im_legend, mylabel, /right, /bottom, box=0, spacing=1.5, $
;           charsize=1.4, psym=psym, color=color, symsize=2.0
       endelse
       if (n_elements(psfile) eq 0) then cc = get_kbrd(1)
    endfor       
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
return
end
