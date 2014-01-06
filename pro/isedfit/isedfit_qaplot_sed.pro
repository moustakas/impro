;+
; NAME:
;   ISEDFIT_QAPLOT_SED
;
; PURPOSE:
;   Generate spectral energy distribution (SED) quality-assurance (QA)
;   plots from the iSEDfit output.
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the QAplots should be
;     written; must match the directory passed to ISEDFIT (default
;     PWD=present working directory) 
;   montegrids_dir - full directory path where the Monte Carlo grids
;     written by ISEDFIT_MONTEGRIDS can be found (default 'montegrids'
;     subdirectory of the PWD=present working directory)
;   index - use this optional input to plot a zero-indexed subset of
;     the full sample (default is to plot everything, although see
;     COMMENTS)
;   galaxy - string array of galaxy names to include in the legend on
;     each page of the QAplot
;   outprefix - optional output prefix string (see ISEDFIT) 
;   pdffile - overwrite the default name of the output PDF file (not
;     typically necessary since the file name matches the ISEDFIT
;     output files); must end with a '.PDF' suffix
;
;   nrandom - build a QAplot for NRANDOM randomly selected galaxies
;     (ignored if INDEX is passed)
;   nsigma - plot photometry detected at less than NSIGMA-sigma as
;     upper limits (default 2.0)
;   xrange, yrange - x- and y-range limits of the plot (useful when
;     you want all the plots to have the same plot limits)
;   xlog - logarithmic wavelength spacing (default linear) 
;
; KEYWORD PARAMETERS:
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;   maxold - see ISEDFIT
;
; OUTPUTS:
;   This routine generates a handy QAplot showing the input
;   photometry, the best-fit (maximum likelihood) SED, and the
;   posterior distributions on several of the key parameters. 
;
; OPTIONAL OUTPUTS:
;   isedfit_results - output data structure containing all the
;     results; see the iSEDfit documentation for a detailed breakdown
;     and explanation of all the outputs  
;
; COMMENTS:
;   This routine should be not be used to plot too many objects,
;   otherwise memory problems may occur; use INDEX or NRANDOM.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Feb 12, U of A
;   jm05aug04uofa - some small changes incorporated
;   jm06mar06uofa - documented; major changes to support the latest
;     version of ISEDFIT 
;   jm06mar23uofa - added FNU, FLAMBDA, MEDIANPLOT, and MAXOLD
;     keywords   
;   jm07jun27nyu  - significantly streamlined; the best-fitting SED
;     models are now read by READ_ISEDFIT()  
;   jm13aug09siena - updated to conform to the latest data model;
;     documentation updated 
;
; Copyright (C) 2005-2007, 2013, John Moustakas
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

pro oplot_posteriors, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte, $
  charsize=charsize, nonorm=nonorm, color_fill=color_fill, $
  color_outline=color_outline, color_monte=color_monte, $
  fill_monte=fill_monte, xlog=xlog, logbins=logbins, _extra=extra, $
  overplot=overplot, nofill=nofill, linestyle=linestyle, thick=thick
; ISEDFIT_QAPLOT_SED: internal support routine
    
    if n_elements(yrange) eq 0 then yrange = [0,1.4]
;   yrange = [0,1.05]

    if (n_elements(xrange) eq 0) then begin
       if keyword_set(overplot) then xrange = !x.crange else $
         xrange = minmax(xx)*[0.95,1.05]
    endif
    if (n_elements(binsize) eq 0) then begin
       if keyword_set(logbins) then $
         binsize = (alog10(xrange[1])-alog10(xrange[0]))/ceil(0.3*sqrt(n_elements(xx))) else $
           binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))
    endif

    if n_elements(color_fill) eq 0 then begin
       if keyword_set(overplot) then color_fill = 'tan' else $
         color_fill = 'powder blue'
    endif
    if n_elements(color_outline) eq 0 then begin
       if keyword_set(overplot) then color_outline = 'tan' else $
         color_outline = 'powder blue'
    endif
    if n_elements(color_monte) eq 0 then color_monte = 'grey80'

    if keyword_set(overplot) eq 0 then begin
;      plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
       plot, [0], [0], xsty=7, ysty=5, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, noerase=noerase, $
         xlog=logbins
    endif
    im_plothist, xx, bin=binsize, peak=1, /overplot, $
      fill=(keyword_set(nofill) eq 0), $
      fcolor=im_color(color_fill), xhist, yhist, charsize=1.0, $
      color=im_color(color_outline,255), logbins=logbins, xlog=logbins, $
      linestyle=linestyle, thick=thick
;   if keyword_set(nomedian) eq 0 then $
;     djs_oplot, median(xx)*[1,1], !y.crange, line=5, thick=8
    if n_elements(monte) ne 0 then begin
       im_plothist, monte, bin=binsize*2, mxhist, myhist, /noplot, $
         _extra=extra, logbins=logbins, xlog=logbins
       if keyword_set(nonorm) then normfactor = 1D else $
         normfactor = max(myhist)/(yrange[1]*0.9)
       if keyword_set(fill_monte) then begin
          im_plothist, monte, bin=binsize*2, /overplot, /fill, $
            normfactor=normfactor, fcolor=im_color(color_monte), $
            color=im_color(color_monte), logbins=logbins, xlog=logbins
       endif else begin
          im_plothist, monte, bin=binsize*2, /overplot, line=1, $
            normfactor=normfactor, color=im_color(color_monte), $
            logbins=logbins, xlog=logbins
;           normfactor=max(myhist)/max(yhist)/1.2
       endelse
    endif
; redraw the axes
    if keyword_set(overplot) then begin
       plot, [0], [0], xsty=3, ysty=1, /noerase, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, ytitle='', $
         xtitle='', xtickinterval=xtickinterval, $
         ytickname=replicate(' ',10), xtickname=replicate(' ',10), $
         charsize=1.0, xlog=logbins, xminor=3, yminor=3, _extra=extra
    endif else begin
;      plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos, $
       plot, [0], [0], xsty=3, ysty=1, /noerase, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, ytitle=ytitle, $
         xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
         ytickname=replicate(' ',10), charsize=1.0, xlog=logbins, $
         xminor=3, yminor=3, _extra=extra
    endelse
return
end

pro isedfit_qaplot_sed, isedfit_paramfile, params=params, thissfhgrid=thissfhgrid, $
  isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, outprefix=outprefix, $
  index=index, galaxy=galaxy1, pdffile=pdffile, nrandom=nrandom, nsigma=nsigma, $
  xrange=in_xrange, yrange=in_yrange, xlog=xlog, isedfit_results=isedfit_results, $
  clobber=clobber

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_qaplot_sed'
       return
    endif

    light = 2.99792458D18       ; speed of light [A/s]
    if n_elements(nsigma) eq 0 then nsigma = 2.0
    
; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()
    if (n_elements(montegrids_dir) eq 0) then montegrids_dir = get_pwd()+'montegrids/'

; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          isedfit_qaplot_sed, params=params[ii], isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, outprefix=outprefix, index=index, $
            galaxy=galaxy1, pdffile=pdffile, xrange=in_xrange, yrange=in_yrange, $
            xlog=xlog, nrandom=nrandom, nsigma=nsigma, clobber=clobber, $
            isedfit_results=isedfit_results1
          isedfit_results1 = struct_trimtags(isedfit_results1,except=['WAVE','FLUX'])
          if ii eq 0 then isedfit_results = isedfit_results1 else $
            isedfit_results = [[isedfit_results],[isedfit_results1]]
       endfor
       return
    endif

; allow the user to overwrite PDFFILE
    fp = isedfit_filepaths(params,outprefix=outprefix,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,sed_pdffile=pdffile)
    if file_test(fp.isedfit_dir+fp.qaplot_sed_pdffile) and $
      keyword_set(clobber) eq 0 then begin
       splog, 'Output file '+fp.qaplot_sed_pdffile+' exists; use /CLOBBER'
       return
    endif

; restore the filters, the best-fitting models, and the posterior
; distributions; optionally restore just NRANDOM random galaxies
    if n_elements(index) eq 0L and n_elements(nrandom) ne 0L then begin
       isedfile = fp.isedfit_dir+fp.isedfit_outfile+'.gz'
       if file_test(isedfile) eq 0 then begin
          splog, 'iSEDfit output file '+isedfile+' not found!'
          return
       endif
       ngal = sxpar(headfits(isedfile,ext=1),'NAXIS2')
       index = shuffle_indx(ngal,num=nrandom)
    endif
    
    isedfit_results = read_isedfit(params=params,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,index=index,outprefix=outprefix,$
      isedfit_post=post,silent=silent,/getmodels)
    ngal = n_elements(isedfit_results)

    filterlist = strtrim(params.filterlist,2)
    nfilt = n_elements(filterlist)
    weff = k_lambda_eff(filterlist=filterlist)
    hwhm = weff*0.0
    
; deal with the galaxy name; note that dimensions of GALAXY must match
; dimensions of INDEX or iSEDFIT
    if n_elements(galaxy1) eq 0L then begin
       mx = max(isedfit_results.isedfit_id)
       fmt = '(I'+string(strlen(strtrim(mx,2)),format='(I0)')+'.'+$
         string(strlen(strtrim(mx,2)),format='(I0)')+')'
;      galaxy = 'Galaxy_'+string(lindgen(ngal),format=fmt)
       galaxy = 'iSEDfit ID '+string(isedfit_results.isedfit_id,format=fmt)
    endif else begin
       if n_elements(index) eq 0L then begin
          if n_elements(galaxy1) ne ngal then begin
             splog, 'Dimensions of GALAXY and ISEDFIT output must agree!'
             return 
          endif
          galaxy = galaxy1
       endif else begin
          galaxy = galaxy1[index]
       endelse
    endelse
    
; generate the QA-plot
    wscale = 1D4
    
    im_plotconfig, 0, pos, psfile=fp.isedfit_dir+fp.qaplot_sed_psfile, $
      ymargin=[1.0,5.9], height=4.1
    pos2 = im_getposition(nx=4,ny=2,ymargin=[6.1,0.8],height=1.8*[1,1],$
      ypage=11.0,yspace=0.5,xmargin=[1.1,0.4],width=1.675*[1,1,1,1],$
      xspace=0.1,xpage=8.5)
    for igal = 0L, ngal-1L do begin
       if ((igal mod 10) eq 0) then print, igal, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'

; redshift and best-fit model
       z = isedfit_results[igal].z
       wave = isedfit_results[igal].wave    ; [A]
       flux = isedfit_results[igal].flux    ; [AB]

; upper panel: SED
       if (isedfit_results[igal].chi2 ge 1E6) then begin
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, yrange=yrange, $
            xtitle='', ytitle='', ytickname=replicate(' ',10), $
            position=pos, xtickname=replicate(' ',10)
          im_legend, ['SED-fitting not possible (see NMINPHOT)'], /left, /top, $
            box=0, spacing=1.5, charsize=1.6
          label = [strtrim(galaxy[igal],2),'z = '+strtrim(string(z,format='(F12.4)'),2)]
          im_legend, label, /right, /bottom, box=0, spacing=1.5, charsize=1.2
       endif else begin
          if (n_elements(in_xrange) eq 2) then xrange1 = in_xrange else begin
             xrange1 = [0.7*min(weff),1.2*max(weff)]
;            xrange1 = [min(filtinfo.weff-1.3*filtinfo.fwhm),$
;              max(filtinfo.weff+2*filtinfo.fwhm)]
             xrange1[0] = xrange1[0]>(700.0*(1+z))
             xrange1[1] = xrange1[1]<max(wave)
             get_element, wave, xrange1, xx
          endelse
          xrange2 = xrange1/(1.0+z)

          notzero = where(isedfit_results[igal].bestmaggies gt 0.0)
          bestmab = -2.5*alog10(isedfit_results[igal].bestmaggies[notzero])

          if (n_elements(in_yrange) eq 2) then yrange = in_yrange else begin
             yrange = fltarr(2)
             yrange[0] = ((max(bestmab)>max(flux[xx[0]:xx[1]]))*1.05)<30.0
             yrange[1] = (min(bestmab)<min(flux[xx[0]:xx[1]]))*0.93
          endelse

          ticks = loglevels(xrange1/wscale)
          djs_plot, [0], [0], /nodata, xrange=xrange1/wscale, yrange=yrange, $
            xsty=9, ysty=1, xtitle='Wavelength \lambda (\mu'+'m)', $
            ytitle='Magnitude (AB)', xlog=xlog, position=pos, $
            xtickv=ticks, xticks=n_elements(ticks)-1

          ticks = loglevels(xrange2/wscale)
          axis, /xaxis, xsty=1, xrange=xrange2/wscale, xlog=xlog, $
            xtitle=textoidl('Rest Wavelength \lambda (\mu'+'m)'), $
            xtickv=ticks, xticks=n_elements(ticks)-1

          djs_oplot, wave/wscale, flux, line=0, color='grey'
          djs_oplot, weff[notzero]/wscale, bestmab, psym=symcat(6,thick=6), $
            symsize=2.5

; overplot the data; distinguish between three different cases, based
; on the input photometry
          mab = maggies2mag(isedfit_results[igal].maggies,nsigma=nsigma,$
            ivar=isedfit_results[igal].ivarmaggies,magerr=maberr,$
            lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
          used = where(mab gt -90.0,nused)
          upper = where(mab lt -90.0 and mabupper gt -90,nupper)

; not used in the fitting          
          notused = where((isedfit_results[igal].maggies gt 0.0) and $ 
            (isedfit_results[igal].ivarmaggies eq 0.0),nnotused)
          nodata = where((isedfit_results[igal].maggies eq 0.0) and $ ; no measurement
            (isedfit_results[igal].ivarmaggies eq 0.0),nnodata)

          if (nused ne 0L) then begin
             oploterror, weff[used]/wscale, mab[used], hwhm[used], $
               mabhierr[used], psym=symcat(16), $
               symsize=1.7, color=im_color('dodger blue'), /hibar, $
               errcolor=im_color('dodger blue'), errthick=!p.thick
             oploterror, weff[used]/wscale, mab[used], hwhm[used], $
               mabloerr[used], psym=3, color=im_color('dodger blue'), /lobar, $
               errcolor=im_color('dodger blue'), errthick=!p.thick
          endif

          if (nupper ne 0) then begin
             djs_oplot, [weff[upper]/wscale], [mabupper[upper]], $
               psym=symcat(11,thick=6), symsize=3.0, color=im_color('forest green')
          endif

          if (nnotused ne 0L) then begin
             mab = maggies2mag(isedfit_results[igal].maggies[notused])
             oploterror, weff[notused]/wscale, mab, hwhm[notused], mab*0.0, $
               psym=symcat(4,thick=6.0), symsize=3.0, color=im_color('firebrick'), $
               errcolor=im_color('firebrick'), errthick=!p.thick
          endif

          if strtrim(params.redcurve,2) eq 'charlot' then begin
             label = [$
               'log (M_{*}/M'+sunsymbol()+') = '+strtrim(string(isedfit_results[igal].mstar,format='(F12.2)'),2),$
               'Age = '+strtrim(string(isedfit_results[igal].age,format='(F12.2)'),2)+' Gyr',$
               '\tau = '+strtrim(string(isedfit_results[igal].tau,format='(F12.1)'),2)+' Gyr',$
               'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit_results[igal].Zmetal/0.019,format='(F12.2)'),2),$
               'A_{V,ISM} = '+strtrim(string(isedfit_results[igal].mu*isedfit_results[igal].av,format='(F12.2)'),2)+' mag',$
               'A_{V,BC} = '+strtrim(string((1-isedfit_results[igal].mu)*isedfit_results[igal].av,format='(F12.2)'),2)+' mag',$
               'log SFR = '+strtrim(string(isedfit_results[igal].sfr,format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}',$
               'log SFR_{100} = '+strtrim(string(isedfit_results[igal].sfr100,format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}',$
               'log b_{100} = '+strtrim(string(isedfit_results[igal].b100,format='(F12.2)'),2)]
;              'log b_{1000} = '+strtrim(string(isedfit_results[igal].b1000,format='(F12.2)'),2)])
          endif else begin
             label = [$
               'log (M_{*}/M'+sunsymbol()+') = '+strtrim(string(isedfit_results[igal].mstar,format='(F12.2)'),2),$
               'Age = '+strtrim(string(isedfit_results[igal].age,format='(F12.2)'),2)+' Gyr',$
               '\tau = '+strtrim(string(isedfit_results[igal].tau,format='(F12.1)'),2)+' Gyr',$
               'Z/Z'+sunsymbol()+' = '+strtrim(string(isedfit_results[igal].Zmetal/0.019,format='(F12.2)'),2),$
               'A_{V} = '+strtrim(string(isedfit_results[igal].av*isedfit_results[igal].mu,format='(F12.2)'),2)+' mag',$
               'log SFR = '+strtrim(string(isedfit_results[igal].sfr,format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}',$
               'log SFR_{100} = '+strtrim(string(isedfit_results[igal].sfr100,format='(F12.2)'),2)+' M'+sunsymbol()+' yr^{-1}',$
               'log b_{100} = '+strtrim(string(isedfit_results[igal].b100,format='(F12.2)'),2)]
;              'log b_{1000} = '+strtrim(string(isedfit_results[igal].b1000,format='(F12.2)'),2)])
          endelse
          
          im_legend, label, /left, /top, box=0, spacing=1.7, charsize=1.2, margin=0
          im_legend, /right, /bottom, box=0, spacing=1.5, charsize=1.2, margin=0, $
            [strtrim(repstr(galaxy[igal],'_',' '),2),$
            'z = '+strtrim(string(z,format='(F12.4)'),2),'\chi_{\nu}^{2} = '+$
            strtrim(string(isedfit_results[igal].chi2,format='(F12.2)'),2)]
       endelse

; lower panels: posterior distributions
       if (isedfit_results[igal].chi2 lt 1E6) then begin

; mass       
          xrange = minmax(post[igal].mstar)*[0.95,1.05]
          oplot_posteriors, post[igal].mstar, pos2[*,0], /noerase, xrange=xrange, $
            xtitle='log (M_{*}/M'+sunsymbol()+')', xtickinterval=0.5

; age       
          xrange = params.age
;         xrange = [min(post[igal].sfrage)<min(post[igal].age),$
;           max(post[igal].sfrage)>max(post[igal].age)]
          oplot_posteriors, post[igal].age, pos2[*,1], /noerase, $
            xtitle='AGE (Gyr)', xtickinterval=1, xrange=xrange
          oplot_posteriors, post[igal].sfrage, pos2[*,1], /overplot
          xyouts, pos2[0,1]+0.01, pos2[3,1]-0.02, 'Age', align=0.0, $
            charsize=1.0, color=im_color('powder blue'), /norm
          xyouts, pos2[0,1]+0.05, pos2[3,1]-0.02, textoidl('<Age>_{SFR}'), $
            align=0.0, charsize=1.0, color=im_color('tan'), /norm

; SFR       
          xrange = [min(post[igal].sfr100)<min(post[igal].sfr),$
            max(post[igal].sfr100)>max(post[igal].sfr)]
;         xrange = [0,(1.1*weighted_quantile(post[igal].sfr100,quant=0.95))>$
;           (1.1*weighted_quantile(post[igal].sfr100,quant=0.95))]
          oplot_posteriors, post[igal].sfr, pos2[*,2], /noerase, xrange=xrange, $
            xtitle='log SFR (M'+sunsymbol()+' yr^{-1})'
          oplot_posteriors, post[igal].sfr100, pos2[*,2], /overplot
          xyouts, pos2[0,2]+0.01, pos2[3,2]-0.02, 'SFR', align=0.0, $
            charsize=1.0, color=im_color('powder blue'), /norm
          xyouts, pos2[0,2]+0.05, pos2[3,2]-0.02, textoidl('<SFR>_{100 Myr}'), $
            align=0.0, charsize=1.0, color=im_color('tan'), /norm

; tau       
          if params.oneovertau then xrange = reverse(1.0/params.tau) else $
            xrange = params.tau
          oplot_posteriors, post[igal].tau, pos2[*,3], /noerase, $
            xrange=xrange, xtitle='\tau (Gyr)', logbins=params.oneovertau

; metallicity       
          oplot_posteriors, post[igal].Zmetal/0.019, pos2[*,4], /noerase, $
            xrange=params.Zmetal/0.019, xtitle='Z/Z'+sunsymbol(), $
            xtickinterval=0.5

; A(V)       
          if strtrim(params.redcurve,2) eq 'none' or max(params.AV) eq 0.0 then begin
             djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,5], $
               xtitle='A_{V} (mag)', charsize=1.0, xtickname=replicate(' ',10), $
               ytickname=replicate(' ',10)
             im_legend, 'Redcurve=none', /left, /top, margin=0, box=0, charsize=1.0
          endif else begin
             if params.flatAV then xrange = params.AV else $
               xrange = [0.0,max(post[igal].AV)*1.1]
             if strtrim(params.redcurve,2) eq 'charlot' then begin
                oplot_posteriors, post[igal].mu*post[igal].AV, pos2[*,5], /noerase, $
                  xrange=xrange, xtitle='A_{V} (mag)'
                oplot_posteriors, (1-post[igal].mu)*post[igal].AV, pos2[*,5], /overplot
                xyouts, pos2[0,5]+0.01, pos2[3,5]-0.02, textoidl('A_{V,ISM}'), align=0.0, $
                  charsize=1.0, color=im_color('powder blue'), /norm
                xyouts, pos2[0,5]+0.06, pos2[3,5]-0.02, textoidl('A_{V,BC}'), $
                  align=0.0, charsize=1.0, color=im_color('tan'), /norm
             endif else begin
                oplot_posteriors, post[igal].AV, pos2[*,5], /noerase, $
                  xrange=xrange, xtitle='A_{V} (mag)'
             endelse
          endelse

; b100 and b1000
          xrange = [min(post[igal].b100)<min(post[igal].b1000),$
            max(post[igal].b100)>max(post[igal].b1000)]
;         xrange = minmax(post[igal].b100)
          oplot_posteriors, post[igal].b100, pos2[*,6], /noerase, xrange=xrange, $
            xtitle='log b'
          oplot_posteriors, post[igal].b1000, pos2[*,6], /overplot
          xyouts, pos2[0,6]+0.01, pos2[3,6]-0.02, textoidl('b_{100}'), align=0.0, $
            charsize=1.0, color=im_color('powder blue'), /norm
          xyouts, pos2[0,6]+0.05, pos2[3,6]-0.02, textoidl('b_{1000}'), $
            align=0.0, charsize=1.0, color=im_color('tan'), /norm
          
; emission lines       
          if params.nebular then begin
             xrange = [0,(1.1*weighted_quantile(post[igal].ewoii,quant=0.95))>$
               (1.1*weighted_quantile(post[igal].ewoiiihb,quant=0.95))>$
               (1.1*weighted_quantile(post[igal].ewniiha,quant=0.95))]
             oplot_posteriors, post[igal].ewoii, pos2[*,7], /noerase, $
               xrange=xrange, xtitle='EW (\AA, rest)'
             oplot_posteriors, post[igal].ewoiiihb, pos2[*,7], /overplot
             oplot_posteriors, post[igal].ewniiha, pos2[*,7], /overplot, color_fill='orange'
             xyouts, pos2[0,7]+0.01, pos2[3,7]-0.02, '[OII]', align=0.0, $
               charsize=1.0, color=im_color('powder blue'), /norm
             xyouts, pos2[0,7]+0.06, pos2[3,7]-0.02, textoidl('[OIII]+H\beta'), $
               align=0.0, charsize=1.0, color=im_color('tan'), /norm
             xyouts, pos2[0,7]+0.01, pos2[3,7]-0.04, textoidl('[NII]+H\alpha'), $
               align=0.0, charsize=1.0, color=im_color('orange'), /norm
          endif else begin
             djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,7], $
               xtitle='EW (\AA, rest)', charsize=1.0, xtickname=replicate(' ',10), $
               ytickname=replicate(' ',10)
             im_legend, 'Nebular=0', /left, /top, margin=0, box=0, charsize=1.0
          endelse

; y-title       
          xyouts, pos2[0,0]-0.05, (pos2[1,0]-pos2[3,4])/2+pos2[3,4], align=0.5, $
            orientation=90, 'Posterior Probability', charsize=1.6, /norm
       endif
    endfor 
    print
    im_plotconfig, psfile=fp.isedfit_dir+fp.qaplot_sed_psfile, /psclose, /pdf
    
return
end
