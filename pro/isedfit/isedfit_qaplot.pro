;+
; NAME:
;   ISEDFIT_QAPLOT
;
; PURPOSE:
;   Generate quality-assurance plots from ISEDFIT output.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   datapath      - should match ISEDFIT [default CWD()] 
;   isedfitprefix - should match ISEDFIT (default 'isedfit')
;
; KEYWORD PARAMETERS:
;   maxold     - see ISEDFIT
;   make_png   - generate PNG output
;   postscript - generate a single postscript output file
;                (stronger than MAKE_PNG)
;
; INPUTS/OUTPUTS:
;   result      - see ISEDFIT
;   result_info - see ISEDFIT
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The cosmological parameters are hard-wired to match
;   ISEDFIT_MODELS.  If POSTSCRIPT=1 then a single, large
;   postscript file is generated for all the objects.  If
;   MAKE_PNG=1 then individual PNG files are generated.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Feb 12, U of A
;   jm05aug04uofa - some small changes incorporated
;   jm06mar06uofa - documented; major changes to support the
;                   latest version of ISEDFIT
;   jm06mar23uofa - added FNU, FLAMBDA, MEDIANPLOT, and MAXOLD
;                   keywords  
;   jm07mar20nyu  - further developments
;   jm07jun27nyu  - significantly streamlined; the best-fitting
;                   SED models are now read by ISEDFIT_RESTORE() 
;
; Copyright (C) 2005-2007, John Moustakas
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

pro isedfit_qaplot, paramfile, isedfit, params=params, super=super, isedpath=isedpath, $
  galaxy=galaxy1, outprefix=outprefix, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
  sfhgrid=sfhgrid, psfile=psfile1, index=index, clobber=clobber, $
  xrange=in_xrange, yrange=in_yrange, xlog=xlog

    light = 2.99792458D18       ; speed of light [A/s]

    nsuper = n_elements(super)
    if nsuper eq 0 or ((n_elements(paramfile) eq 0) and $
      (n_elements(params) eq 0)) then begin
       doc_library, 'isedfit_qaplot'
       return
    endif

; read the parameter file and parse to get the relevant path and
; filenames; optionally overwrite SFHGRID in the parameter file
    if (n_elements(isedpath) eq 0) then isedpath = './'
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(paramfile,sfhgrid=sfhgrid)
    
; SUPER can be a vector
    if nsuper gt 1 then begin
       for ii = 0, nsuper-1 do begin
          isedfit_qaplot, params=params, super=super, isedpath=isedpath, galaxy=galaxy1, $
            outprefix=outprefix, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            psfile=psfile1, index=index, clobber=clobber, xrange=xrange1, $
            yrange=yrange, xlog=xlog
       endfor
       return
    endif

; allow the user to overwrite PSFILE
    fp = isedfit_filepaths(params,super=super,outprefix=outprefix,$
      isedpath=isedpath,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)
    if (n_elements(psfile1) eq 0) then $
      psfile = isedpath+strtrim(fp.qaplot_psfile,2) else $
        psfile = psfile1
    if file_test(psfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+psfile+' exists; use /CLOBBER'
       return
    endif

; restore the filters and best-fitting models
    filterlist = strtrim(params.filterlist,2)
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilt = n_elements(filterlist)

    model = isedfit_restore(paramfile,isedfit,params=params,super=super,$
      isedpath=isedpath,index=index,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,$
      outprefix=outprefix,silent=silent)
    ngal = n_elements(isedfit)

; deal with the galaxy name; note that dimensions of GALAXY must match
; dimensions of INDEX or iSEDFIT
    if (n_elements(galaxy1) eq 0L) then begin
       fmt = '(I'+string(5L,format='(I0)')+'.'+string(5L,format='(I0)')+')'
       galaxy = 'Galaxy_'+string(lindgen(ngal),format=fmt)
    endif else begin
       if (n_elements(index) eq 0L) then begin
          if (n_elements(galaxy1) ne ngal) then begin
             splog, 'GALAXY and ISEDFIT output must have '+$
               'the same number of elements'
             return
          endif
          galaxy = galaxy1
       endif else begin
          if (n_elements(galaxy1) ne n_elements(index)) then begin
             splog, 'GALAXY and INDEX  must have '+$
               'the same number of elements'
             return
          endif
       endelse
       galaxy = galaxy1
    endelse 
    
; generate the QA-plot
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    im_plotconfig, 0, pos, psfile=psfile, ymargin=[1.0,1.1], height=5
;   for igal = 54, ngal-1L do begin
    for igal = 0L, ngal-1L do begin
       if ((igal mod 10) eq 0) then print, igal, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'

; redshift and best-fit model
       zobj = isedfit[igal].zobj
       wave = model[igal].wave    ; [A]
       flux = model[igal].flux    ; [AB]

; make the plot       
       if (isedfit[igal].chi2 ge 1E6) then begin
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, yrange=yrange, $
            xtitle=xtitle1, ytitle=ytitle1, ytickname=replicate(' ',10), $
            position=pos        ; ymargin=[4,3], 
          legend, ['No mass estimate available'], /left, /top, $
            box=0, spacing=1.5, charsize=1.6
          label = [strtrim(galaxy[igal],2),'z = '+strtrim(string(zobj,format='(F12.4)'),2)]
          legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.6
       endif else begin
          if (n_elements(in_xrange) eq 2) then xrange1 = in_xrange else begin
             xrange1 = [min(filtinfo.weff-1.3*filtinfo.fwhm),$
               max(filtinfo.weff+2*filtinfo.fwhm)]
             xrange1[0] = xrange1[0]>(700.0*(1+zobj))
             xrange1[1] = xrange1[1]<max(wave)
             get_element, wave, xrange1, xx
          endelse
          xrange2 = xrange1/(1.0+zobj)

          weff = filtinfo.weff
          hwhm = filtinfo.fwhm/2.0

          notzero = where(isedfit[igal].bestmaggies gt 0.0)
          bestmab = -2.5*alog10(isedfit[igal].bestmaggies[notzero])

          if (n_elements(in_yrange) eq 2) then yrange = in_yrange else begin
             yrange = fltarr(2)
             yrange[0] = ((max(bestmab)>max(flux[xx[0]:xx[1]]))*1.05)<30.0
             yrange[1] = (min(bestmab)<min(flux[xx[0]:xx[1]]))*0.93
          endelse

          djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
            xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=xlog, $
            xtickinterval=2000, position=pos, xtickformat='(I0)'
          axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, $
            xlog=xlog

          djs_oplot, wave, flux, line=0, color='grey'
          djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), symsize=2.5

; overplot the data; distinguish between three different cases, based
; on the input photometry
          mab = maggies2mag(isedfit[igal].maggies,ivar=isedfit[igal].ivarmaggies,$
            magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper,nsigma=2.0)
          used = where(mab gt -90.0,nused)
          upper = where(mab lt -90.0 and mabupper gt -90,nupper)

;         used = where((isedfit[igal].maggies gt 0.0) and $ ; used in the fitting
;           (isedfit[igal].ivarmaggies gt 0.0),nused)
;         upper = where((isedfit[igal].maggies le 0.0) and $ ; upper limit
;           (isedfit[igal].ivarmaggies gt 0.0),nupper)
          notused = where((isedfit[igal].maggies gt 0.0) and $ ; not used in the fitting
            (isedfit[igal].ivarmaggies eq 0.0),nnotused)
          nodata = where((isedfit[igal].maggies eq 0.0) and $ ; no measurement
            (isedfit[igal].ivarmaggies eq 0.0),nnodata)

          if (nused ne 0L) then begin
             oploterror, weff[used], mab[used], hwhm[used], mabhierr[used], psym=symcat(16), $
               symsize=2.0, color=im_color('firebrick'), /hibar, $
               errcolor=im_color('firebrick'), errthick=!p.thick
             oploterror, weff[used], mab[used], hwhm[used], mabloerr[used], psym=3, $
               color=im_color('firebrick'), /lobar, $
               errcolor=im_color('firebrick'), errthick=!p.thick
          endif

          if (nupper ne 0) then begin
             djs_oplot, weff[upper], mabupper[upper], $
               psym=symcat(11,thick=4), symsize=3.0, color=im_color('dodger blue')
;            oploterror, weff[upper], mabupper[upper], hwhm[upper], mabupper[upper]*0, $
;              psym=symcat(11,thick=4), symsize=3.0, color=im_color('blue'), $
;              errcolor=im_color('blue'), errthick=!p.thick
          endif

          if (nnotused ne 0L) then begin
             mab = maggies2mag(isedfit[igal].maggies[notused])
             oploterror, weff[notused], mab, hwhm[notused], mab*0.0, $
               psym=symcat(4,thick=6.0), symsize=3.0, color=im_color('forest green'), $
               errcolor=im_color('forest green'), errthick=!p.thick
          endif

; legend
          label = [$
            'log (M/M'+sunsymbol()+') = '+strtrim(string(isedfit[igal].mass,format='(F12.2)'),2)+$
            ' ('+strtrim(string(isedfit[igal].mass_50,format='(F12.2)'),2)+')',$
            'A_{V} = '+strtrim(string(isedfit[igal].av*isedfit[igal].mu,format='(F12.3)'),2)+$
            ' ('+strtrim(string(isedfit[igal].av_50*isedfit[igal].mu_50,format='(F12.3)'),2)+')',$
            'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit[igal].Z/0.02,format='(F12.2)'),2)+$
            ' ('+strtrim(string(isedfit[igal].Z_50/0.02,format='(F12.2)'),2)+')',$
            '\tau = '+strtrim(string(isedfit[igal].tau,format='(F12.1)'),2)+$
            ' ('+strtrim(string(isedfit[igal].tau_50,format='(F12.2)'),2)+') Gyr',$
            'Age = '+strtrim(string(isedfit[igal].age,format='(F12.2)'),2)+$
            ' ('+strtrim(string(isedfit[igal].age_50,format='(F12.2)'),2)+') Gyr',$
            'log \psi_{100} = '+strtrim(string(isedfit[igal].sfr100,format='(F12.3)'),2)+$
            ' ('+strtrim(string(isedfit[igal].sfr100_50,format='(F12.3)'),2)+') M'+sunsymbol()+' yr^{-1}',$
            'b_{100} = '+strtrim(string(isedfit[igal].b100,format='(F12.3)'),2)+$
            ' ('+strtrim(string(isedfit[igal].b100_50,format='(F12.3)'),2)+')']
          legend, textoidl(label), /left, /top, box=0, spacing=1.7, charsize=1.1, margin=0
          label = textoidl([strtrim(repstr(galaxy[igal],'_',' '),2),$
            'z = '+strtrim(string(zobj,format='(F12.4)'),2),'\chi_{\nu}^{2} = '+$
            strtrim(string(isedfit[igal].chi2,format='(F12.2)'),2)])
          legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.4, margin=0
       endelse
    endfor 
    im_plotconfig, psfile=psfile, /psclose, /pdf;, /gzip
    
return
end
