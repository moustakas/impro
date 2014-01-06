;+
; NAME:
;   ISEDFIT_MONTEGRIDS
;
; PURPOSE:
;   Build the Monte Carlo star formation history grids.
;
; INPUTS: 
;   isedfit_paramfile - iSEDfit parameter file
;
; OPTIONAL INPUTS: 
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the iSEDfit models and
;     output files will be written (default PWD=present working
;     directory) 
;   montegrids_dir - full directory path where the Monte Carlo grids
;     should be written (default 'montegrids' subdirectory of the
;     PWD=present working directory) 
;
;   modelchunksize, minichunksize - to deal with disk space issues
;     this routine splits the full set of models (NMODEL) into
;     MODELCHUNKSIZE sized chunks (as specified in the
;     ISEDFIT_PARAMFILE parameter file); furthemore, to circumvent
;     additional possible memory issues, each chunk is further divided
;     into MINICHUNKSIZE (default 500) sized chunks; all these tricks 
;     *should be* transparent to the user and are only documented here
;     for completeness 
;
; KEYWORD PARAMETERS: 
;   clobber - delete old files from a previous call to this routine,
;     including all the Monte Carlo distributions of parameter values
;     (the user is prompted to confirm the clobber to avoid
;     inadvertently deleting a large grid!)
;   debug - make some debugging plots and wait for a keystroke
;     (currently not a very useful switch)
;
; OUTPUTS: 
;   The grids (spectra, parameters, etc.) are written out in a data
;   model that iSEDfit understands, and which should be transparent to
;   the user.  However, this routine also generates a QAplot (written
;   to the ISEDFIT_DIR directory) which should be inspected to ensure
;   that proper parameter priors have been chosen.  
;
; COMMENTS:
;   Some ToDo items:
;     * output the UV slope (beta) and maybe D(4000)
;     * better debugging diagnostic plots 
;     * Speed the routine up?
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 20, NYU - written
;   jm09may22nyu - generalized how the reddening curves are treated 
;   jm10jul11ucsd - include SFR measures on various timescales and the
;     birthrate parameter at each age
;   jm10nov12ucsd - streamlined and parameter file adopted
;   jm11aug24ucsd - documentation updated
;   jm13jan13siena - semi-major changes in preparation for public
;     release 
;   jm13aug09siena - significant update and rewrite to conform to a
;     new and much simpler data model
;
; Copyright (C) 2009-2011, 2013, John Moustakas
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

function init_montegrid, nmodel, nmaxburst=nmaxburst
; internal support routine: initialize the structure containing the
; parameters of the Monte Carlo grid 

    burstarray1 = -1.0
    if (nmaxburst gt 1) then burstarray1 = fltarr(nmaxburst)-1.0
    
    montegrid = {$
      chunkindx:                 -1,$
      modelindx:                 -1,$
      age:                     -1.0,$
      tau:                     -1.0,$
      Zmetal:                  -1.0,$
      AV:                       0.0,$ ; always initialize with zero to accommodate the dust-free models!
      mu:                       1.0,$ ; always default to 1.0!
      oiiihb:                  -1.0,$ ; [OIII] 5007/H-beta
      nlyc:                    -1.0,$ ; number of Lyman-continuum photons
      mstar:                   -1.0,$
      sfr:                     -1.0,$
      sfr100:                  -1.0,$
      b100:                    -1.0,$
      b1000:                   -1.0,$
      sfrage:                  -1.0,$
      nburst:                     0,$
      trunctau:                -1.0,$ ; burst truncation time scale
      tburst:           burstarray1,$
      dtburst:          burstarray1,$ ; 
      fburst:           burstarray1}  ; burst mass fraction
    montegrid = replicate(montegrid,nmodel)

return, montegrid
end

function get_bracket_files, info, sspinfo=sspinfo
; internal support routine: get the two SSPs that bracket the desired
; metallicity, to allow for interpolation
    ninfo = n_elements(info)
    bracket_files = strarr(2,ninfo)

    for ii = 0L, ninfo-1 do begin
       indx = findex(sspinfo.Zmetal,info[ii].Zmetal)
       below = fix(indx)
       above = ceil(indx)
       bracket_files[0,ii] = sspinfo.sspfile[below]
       bracket_files[1,ii] = sspinfo.sspfile[above]
    endfor
    
return, bracket_files
end

function read_and_interpolate, files, info=info, $
  ssppath=ssppath, fits_grid=fits_grid
; internal support routine: read and interpolate the base model file
; onto the desired metallicity grid 

    files = strtrim(files,2)
    junk1 = file_search(ssppath+files[0],count=c1)
    junk2 = file_search(ssppath+files[1],count=c2)
    if (c1 eq 0) or (c2 eq 0) then message, 'Problem finding files'
    
    nobj = n_elements(info)
    if (strtrim(files[0],2) eq strtrim(files[1],2)) then begin
       fits = gz_mrdfits(ssppath+files[0],1,/silent)
       fits = replicate(temporary(fits),nobj) ; identical for all objects
    endif else begin
       fits_grid = [gz_mrdfits(ssppath+files[0],1,/silent),$
         gz_mrdfits(ssppath+files[1],1,/silent)]
       fits = replicate(fits_grid[0],nobj)
       indx = findex(fits_grid.Zmetal,info.Zmetal)
; this interpolate can be very slow because we're building an
; NPIX,NAGE,NOBJ elements array
;      t0 = systime(1)
       fits.flux = interpolate(fits_grid.flux,indx)
;      splog, 'Total time (sec) = ', (systime(1)-t0)
       fits.mstar = interpolate(fits_grid.mstar,indx)
       fits.Zmetal = interpolate(fits_grid.Zmetal,indx)
    endelse       

return, fits
end

pro oplot_priors, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, charsize=charsize, nonorm=nonorm, $
  color_fill=color_fill, color_outline=color_outline, $
  xlog=xlog, logbins=logbins, _extra=extra, edge=edge, $
  overplot=overplot, nofill=nofill, linestyle=linestyle, thick=thick
; ISEDFIT_QAPLOT_SED: internal support routine
    
    if n_elements(yrange) eq 0 then yrange = [0,1.4]
    if (n_elements(xrange) eq 0) then begin
       if keyword_set(overplot) then xrange = !x.crange else $
         xrange = minmax(xx)*[0.95,1.05]
    endif
    if (n_elements(binsize) eq 0) then begin
       if keyword_set(logbins) then $
         binsize = (alog10(xrange[1])-alog10(xrange[0]))/$
         ceil(0.3*sqrt(n_elements(xx))) else $
           binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))
    endif

    csize = 1.4
    
    if n_elements(color_fill) eq 0 then begin
       if keyword_set(overplot) then color_fill = 'tan' else $
         color_fill = 'powder blue'
    endif
    if n_elements(color_outline) eq 0 then begin
       if keyword_set(overplot) then color_outline = 'tan' else $
         color_outline = 'powder blue'
    endif

    if keyword_set(overplot) eq 0 then begin
       plot, [0], [0], xsty=7, ysty=5, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, noerase=noerase, $
         xlog=logbins
    endif
    im_plothist, xx, bin=binsize, peak=1, /overplot, $
      fill=(keyword_set(nofill) eq 0), edge=edge, $
      fcolor=im_color(color_fill), xhist, yhist, charsize=csize, $
      color=im_color(color_outline,255), logbins=logbins, xlog=logbins, $
      linestyle=linestyle, thick=thick
; redraw the axes
    if keyword_set(overplot) then begin
       plot, [0], [0], xsty=3, ysty=1, /noerase, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, ytitle='', $
         xtitle='', xtickinterval=xtickinterval, $
         ytickname=replicate(' ',10), xtickname=replicate(' ',10), $
         charsize=csize, xlog=logbins, xminor=3, yminor=3, _extra=extra
    endif else begin
       plot, [0], [0], xsty=3, ysty=1, /noerase, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, ytitle=ytitle, $
         xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
         ytickname=replicate(' ',10), charsize=csize, xlog=logbins, $
         xminor=3, yminor=3, _extra=extra
    endelse
return
end

pro montegrids_qaplot, montegrid, params=params, qafile=qafile
; internal support routine: build a QAplot 
    
; page 1: parameter choices
    im_plotconfig, 0, pos, psfile=qafile, ymargin=[0.4,1.0], height=7.8, $
      xmargin=[0.75,0.75], width=7.0

    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5
    label = ['iSEDfit Prior Parameters',' ',$
      'prefix='+params.prefix,$
      'spsmodels='+params.spsmodels,$
      'imf='+params.imf,$
      'redcurve='+params.redcurve,$
      ' ',$
      'h100='+string(params.h100,format='(G0)'),$
      'omega0='+string(params.omega0,format='(G0)'),$
      'omegal='+string(params.omegal,format='(G0)'),$
      ' ',$
      'nmodel='+string(params.nmodel,format='(G0)'),$
      'ndraw='+string(params.ndraw,format='(G0)'),$
      'nminphot='+string(params.nminphot,format='(G0)'),$
      'galchunksize='+string(params.galchunksize,format='(G0)'),$
      'modelchunksize='+string(params.modelchunksize,format='(G0)'),$
      ' ',$
      'pburst='+string(params.pburst,format='(G0)'),$
      'interval_pburst='+strtrim(string(params.interval_pburst,format='(F12.2)'),2)+' Gyr',$
      'nmaxburst='+string(params.nmaxburst,format='(G0)'),$
      'fractrunc='+string(params.fractrunc,format='(G0)'),$
      ' ',$
      'igm='+string(params.igm,format='(G0)'),$
      'nebular='+string(params.nebular,format='(G0)'),$
      'oneovertau='+string(params.oneovertau,format='(G0)'),$
      'delayed='+string(params.delayed,format='(G0)'),$
      'bursttype='+string(params.bursttype,format='(G0)'),$
      ' ',$
      'flatav='+string(params.flatav,format='(G0)'),$
      'flatmu='+string(params.flatmu,format='(G0)'),$
      'flatfburst='+string(params.flatfburst,format='(G0)'),$
      'flatdtburst='+string(params.flatdtburst,format='(G0)'),$
      ' ',$
      'user_redshift='+string(params.user_redshift,format='(G0)'),$
      'zminmax='+'['+strjoin(string(params.zminmax,format='(G0.0)'),',')+']',$
      'zbin='+string(params.zbin,format='(G0)'),$
      'nzz='+string(params.nzz,format='(G0)'),$
      'zlog='+string(params.zlog,format='(I0)')]
    im_legend, label, /left, /top, box=0, charsize=1.3, margin=0, /notextoidl

    csize = 1.4
    lsize = 1.2
    
; page 2 - primary priors
    pos = im_getposition(nx=2,ny=3,xmargin=[0.75,0.75],ymargin=[0.4,1.0],$
      width=3.2*[1,1],height=2.6*[1,1,1],xspace=0.6,yspace=0.9*[1,1],$
      xpage=8.5,ypage=11.0,landscape=landscape)

; age
    xrange = minmax(montegrid.age)
    oplot_priors, montegrid.age, pos[*,0], xrange=xrange, $
      xtitle='Age (Gyr)'
    im_legend, '['+strjoin(string(params.age,format='(G0.0)'),',')+'] Gyr', $
      /left, /top, box=0, charsize=lsize, margin=0

    xyouts, pos[0,0], pos[3,0]+0.04, 'Input (Primary) Priors', $
      align=0.0, /normal, charsize=1.5
    
; tau
    if params.oneovertau then begin
       xrange = reverse(1.0/params.tau)
       oplot_priors, montegrid.tau, pos[*,1], /noerase, xrange=xrange, $
         xtitle='\tau (Gyr)', /logbins
;      oplot_priors, 1.0/montegrid.tau, pos[*,1], /noerase, xrange=xrange, $
;        xtitle='1/\tau (Gyr^{-1})';, /logbins
       im_legend, ['Oneovertau=1','['+strjoin(string(params.tau,format='(G0.0)'),',')+'] Gyr^{-1}'], $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       xrange = params.tau
       oplot_priors, montegrid.tau, pos[*,1], /noerase, xrange=xrange, $
         xtitle='\tau (Gyr)'
       im_legend, '['+strjoin(string(params.tau,format='(G0.0)'),',')+'] Gyr', $
         /left, /top, box=0, charsize=lsize, margin=0
    endelse

; metallicity       
    oplot_priors, montegrid.Zmetal, pos[*,2], /noerase, $
      xrange=params.Zmetal, xtitle='Z_{metal}' ; /Z'+sunsymbol()
    im_legend, '['+strjoin(string(params.Zmetal,format='(G0.0)'),',')+']', $
      /left, /top, box=0, charsize=lsize, margin=0
;   im_legend, 'Z'+sunsymbol()+'=0.019', /left, /top, box=0, charsize=csize, margin=0

; A(V)       
    if strtrim(params.redcurve,2) eq 'none' then begin
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,3], $
         xtitle='A_{V} (mag)', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'Redcurve=none', /left, /top, margin=0, box=0, charsize=lsize
    endif else begin
       if params.flatAV then xrange = params.AV else $
         xrange = [0.0,max(montegrid.AV)*1.1]
       oplot_priors, montegrid.AV, pos[*,3], /noerase, $
         xrange=xrange, xtitle='A_{V} (mag)'
       if params.flatAV then begin
          im_legend, ['Redcurve='+strtrim(params.redcurve,2),$
            '['+strjoin(string(params.AV,format='(G0.0)'),',')+']'], $
            /left, /top, box=0, charsize=lsize, margin=0
       endif else begin
          im_legend, ['Redcurve='+strtrim(params.redcurve,2),$
            '['+string(params.AV[0],format='(G0.0)')+'*\Gamma('+$
            string(params.AV[1],format='(G0.0)')+')] mag'], $
            /left, /top, box=0, charsize=lsize, margin=0
       endelse
    endelse

; mu
    if strtrim(params.redcurve,2) eq 'charlot' then begin
       if params.flatmu then xrange = params.mu else $
         xrange = [0.0,max(montegrid.mu)*1.1]
       oplot_priors, montegrid.mu, pos[*,4], /noerase, $
         xrange=xrange, xtitle='\mu'
       if params.flatmu then begin
          im_legend, ['Redcurve='+strtrim(params.redcurve,2),$
            '['+strjoin(string(params.mu,format='(G0.0)'),',')+']'], $
            /left, /top, box=0, charsize=lsize, margin=0
       endif else begin
          im_legend, ['Redcurve='+strtrim(params.redcurve,2),$
            '['+string(params.mu[0],format='(G0.0)')+'*\Gamma('+$
            string(params.mu[1],format='(G0.0)')+')]'], $
            /left, /top, box=0, charsize=lsize, margin=0
       endelse
    endif else begin
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,4], $
         xtitle='\mu', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'Redcurve='+strtrim(params.redcurve,2), $
         /left, /top, margin=0, box=0, charsize=lsize
    endelse

; [OIII]/H-beta
    if params.nebular then begin
       oplot_priors, montegrid.oiiihb, pos[*,5], /noerase, $
         xrange=params.oiiihb, xtitle='log_{10}([OIII \lambda5007]/H\beta)'
       im_legend, '['+strjoin(string(params.oiiihb,format='(G0.0)'),',')+'] dex', $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,5], $
         xtitle='log_{10} ([OIII \lambda5007]/H\beta)', $
         charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'Nebular=0', /left, /top, margin=0, box=0, charsize=lsize
    endelse

; page 3: burst/truncation parameters

; nburst
    if params.pburst gt 0.0 then begin
       oplot_priors, montegrid.nburst, pos[*,0], xtitle='N_{burst}', $
         binsize=0.5, edge=-1, xtickinterval=1
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,0], $
         xtitle='N_{burst}', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'No bursts', /left, /top, margin=0, box=0, charsize=lsize
    endelse

    xyouts, pos[0,0], pos[3,0]+0.04, 'Input Burst/Truncated SFH Priors', $
      align=0.0, /normal, charsize=1.5

    xyouts, pos[0,0], pos[3,0]+0.04, 'Input Burst/Truncated SFH Priors', $
      align=0.0, /normal, charsize=1.5

; tburst
    if params.pburst gt 0.0 then begin
       yes = where(montegrid.tburst gt 0.0)
       xrange = minmax((montegrid.tburst)[yes])
       oplot_priors, (montegrid.tburst)[yes], pos[*,1], /noerase, xrange=xrange, $
         xtitle='t_{burst} (Gyr)'
       im_legend, '['+strjoin(string(params.tburst,format='(G0.0)'),',')+'] Gyr', $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, $
         xtitle='t_{burst} (Gyr)', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'No bursts', /left, /top, margin=0, box=0, charsize=lsize
    endelse

; fburst
    if params.pburst gt 0.0 then begin
       yes = where(montegrid.fburst gt 0.0)
       xrange = minmax((montegrid.fburst)[yes])
       oplot_priors, (montegrid.fburst)[yes], pos[*,2], /noerase, xrange=xrange, $
         xtitle='f_{burst}';, logbins=params.flatfburst eq 0
       im_legend, ['Flatfburst='+strtrim(params.flatfburst,2),$
         '['+strjoin(string(params.fburst,format='(G0.0)'),',')+']'], $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,2], /noerase, $
         xtitle='f_{burst}', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'No bursts', /left, /top, margin=0, box=0, charsize=lsize
    endelse

; dtburst
    if params.pburst gt 0.0 then begin
       yes = where(montegrid.dtburst gt 0.0)
       xrange = minmax((montegrid.dtburst)[yes])
       oplot_priors, (montegrid.dtburst)[yes], pos[*,3], /noerase, xrange=xrange, $
         xtitle='\Delta'+'t_{burst} (Gyr)';, logbins=params.flatdtburst eq 0
       im_legend, ['Flatdtburst='+strtrim(params.flatdtburst,2),$
         '['+strjoin(string(params.dtburst,format='(G0.0)'),',')+'] Gyr'], $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,3], /noerase, $
         xtitle='\Delta'+'t_{burst} (Gyr)', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'No bursts', /left, /top, margin=0, box=0, charsize=lsize
    endelse

; trunctau
    if params.fractrunc gt 0.0 then begin
       yes = where(montegrid.trunctau gt 0.0)
       xrange = minmax(montegrid[yes].trunctau)
       oplot_priors, montegrid[yes].trunctau, pos[*,4], /noerase, xrange=xrange, $
         xtitle='\tau_{trunc} (Gyr)'
       im_legend, ['Fractrunc='+string(params.fractrunc,format='(G0.0)'),$
         '['+strjoin(string(params.trunctau,format='(G0.0)'),',')+'] Gyr'], $
         /left, /top, box=0, charsize=lsize, margin=0
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,4], /noerase, $
         xtitle='\tau_{trunc} (Gyr)', charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'No truncated SFHs', /left, /top, margin=0, box=0, charsize=lsize
    endelse

; page 4: derived parameters

; t/tau
    good = where(montegrid.tau gt 0.0,ngood)
    if ngood gt 0L then begin
       xrange = minmax(montegrid[good].age/montegrid[good].tau)
       oplot_priors, montegrid[good].age/montegrid[good].tau, pos[*,0], $
         xrange=xrange, xtitle='Age/\tau', /logbins
    endif else begin
       djs_plot, [0], [0], /nodata, position=pos[*,0], $
         xtitle='Age/\tau', $
         charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'tau=0', /left, /top, margin=0, box=0, charsize=lsize
    endelse

    xyouts, pos[0,0], pos[3,0]+0.04, 'Derived (Secondary) Priors', $
      align=0.0, /normal, charsize=1.5
    
; SFR-weighted age    
    xrange = minmax(montegrid.sfrage)
    oplot_priors, montegrid.sfrage, pos[*,1], /noerase, $
      xrange=xrange, xtitle='Age (SFR-weighted, Gyr)'

; b100
    xrange = minmax(montegrid.b100)
    oplot_priors, montegrid.b100, pos[*,2], /noerase, xrange=xrange, $
      xtitle='log b_{100}'

; b1000
    xrange = minmax(montegrid.b1000)
    oplot_priors, montegrid.b1000, pos[*,3], /noerase, xrange=xrange, $
      xtitle='log b_{1000}'

; EW([OIII]+Hbeta)
    if params.nebular then begin
       xrange = [0,(1.1*weighted_quantile(montegrid.ewoii,quant=0.95))>$
         (1.1*weighted_quantile(montegrid.ewoiiihb,quant=0.95))>$
         (1.1*weighted_quantile(montegrid.ewniiha,quant=0.95))]
       oplot_priors, montegrid.ewoiiihb, pos[*,4], /noerase, $
         xtitle='EW (\AA, rest)', logbins=0, $
         /nofill, line=0, xrange=xrange, thick=8
       oplot_priors, montegrid.ewniiha, pos[*,4], /overplot, logbins=0, $
         /nofill, line=5, color_outline='tan', xrange=xrange, thick=8
       oplot_priors, montegrid.ewoii, pos[*,4], /overplot, logbins=0, $
         /nofill, line=3, color_outline='orange', xrange=xrange, thick=8
       im_legend, ['EW([OIII]+H\beta)','EW([NII]+H\alpha)','EW([OII])'], $
         /right, /top, box=0, charsize=1.1, line=[0,5,3], pspacing=1.9, $
         margin=0, color=['powder blue','tan','orange'], $
         textcolor=['powder blue','tan','orange']
    endif else begin
       djs_plot, [0], [0], /nodata, /noerase, position=pos[*,4], $
         xtitle='EW (\AA, rest)', $
         charsize=csize, xtickname=replicate(' ',10), $
         ytickname=replicate(' ',10)
       im_legend, 'Nebular=0', /left, /top, margin=0, box=0, charsize=lsize       
    endelse
    
;; SFR/SFR100
;    mm = 1D-14
;    xrange = minmax((montegrid.sfr>mm)*1D11)
;    oplot_priors, (montegrid.sfr>mm)*1D11, pos[*,2], /noerase, xrange=xrange, $
;      xtitle='SFR (10^{-11} M'+sunsymbol()+' yr^{-1})', /logbins
;;   oplot_priors, montegrid.sfr100*1D11, pos[*,2], /overplot

    im_plotconfig, psfile=qafile, /psclose, /pdf

return    
end   

function build_modelgrid, montegrid, params=params, debug=debug, $
  sspinfo=sspinfo, ssppath=ssppath, minichunksize=minichunksize, $
  ichunk=ichunk, nchunk=nchunk
; internal support routine: workhorse code which builds the spectra  

    if (n_elements(tbc) eq 0) then tbc = 0.01D9 ; dust dispersal timescale [10 Myr]

; define some reddening curve variables
    case strlowcase(strtrim(params.redcurve,2)) of
       'none': 
       'calzetti': calzetti = 1
       'charlot': charlot = 1
       'odonnell': odonnell = 1
       'smc': smc = 1
       else: message, 'Unrecognized reddening curve'
    endcase

; loop on each mini-chunk of models
    nmodel = n_elements(montegrid)
    nmini = ceil(nmodel/float(minichunksize))
    
    for imini = 0, nmini-1 do begin
       i1 = imini*minichunksize
       i2 = ((imini*minichunksize+minichunksize)<nmodel)-1L
       nthese = i2-i1+1L & these = lindgen(nthese)+i1

       modelgrid1 = montegrid[these]
       
; get the two SSPs that bracket the desired metallicity, to allow for
; interpolation; do this here so that we're not reading in the same
; set of models multiple times; jm13aug07siena: interpolating a large
; number of files is very slow and memory-intensive, so subdivide the
; SSPs into mini-chunks
       bracket_sspfile = get_bracket_files(modelgrid1,sspinfo=sspinfo)
       uindx = uniq(bracket_sspfile[0,*],sort(bracket_sspfile[0,*]))
       ubracket_sspfile = bracket_sspfile[*,uindx]
       nssp = n_elements(uindx)

       for issp = 0, nssp-1 do begin ; loop on SSP metallicity
          sspindx = where(ubracket_sspfile[0,issp] eq bracket_sspfile[0,*],nsspindx)
          modelgrid1[sspindx].modelindx = sspindx+i1

          sspfits = read_and_interpolate(ubracket_sspfile[*,issp],$
            info=modelgrid1[sspindx],ssppath=ssppath)

; build an emission-line friendly wavelength array; this should be
; identical for all SSPs
          if params.nebular then begin
             junk = isedfit_nebular(1D,wave=emwave,vsigma=vsigma,$
               inst_vsigma=sspinfo.inst_vsigma)

             wave = [emwave,sspfits[0].wave]
             wave = float(wave[sort(wave)])
             if nsspindx eq 1L then begin
                sspflux = interpolate(sspfits.flux,findex(sspfits[0].wave,wave),$
                  findgen(n_elements(sspfits[0].age)),/grid)
             endif else begin
                sspflux = interpolate(sspfits.flux,findex(sspfits[0].wave,wave),$
                  findgen(n_elements(sspfits[0].age)),findgen(nsspindx),/grid)
             endelse
          endif else begin
             wave = sspfits[0].wave
             sspflux = sspfits.flux
          endelse
          npix = n_elements(wave)

; initialize the output SED data structure; we waste a little disk
; space by replicating the wavelength array for each model, but
; it's fine, disk is cheap
          if issp eq 0 then begin
             spectags = replicate({ewoii: -1.0, ewoiiihb: -1.0, ewniiha: -1.0, $
               wave: float(wave), flux: fltarr(npix)},nthese)
             modelgrid1 = struct_addtags(temporary(modelgrid1),spectags)
          endif

          delvarx, rv
          klam = k_lambda(wave,r_v=rv,charlot=charlot,calzetti=calzetti,$
            odonnell=odonnell,smc=smc,/silent)
          klam = rebin(reform(klam,npix,1),npix,n_elements(sspfits[0].age)) ; [NPIX,NAGE]

; convolve each SSP with the specified SFH
          for jj = 0, nsspindx-1 do begin
             if ((jj mod 10) eq 0) then print, format='("Chunk=",I4.4,"/",I4.4,", '+$
               'Minichunk=",I4.4,"/",I4.4,", SSP(Z)=",I3.3,"/",I3.3,", '+$
               'Model=",I4.4,"/",I4.4,"   ",A5,$)', ichunk+1, nchunk, $
               imini+1, nmini, issp+1, nssp, jj, nsspindx, string(13b)
; attenuation; for the Charlot & Fall attenuation curve adopt the
; time-dependent dust model of Pacifici et al. 2012
             if keyword_set(charlot) then begin ; special case
                alam_bc = (1.0-modelgrid1[sspindx[jj]].mu)*modelgrid1[sspindx[jj]].av*$
                  (wave/5500D)^(-1.3)
                alam_ism = modelgrid1[sspindx[jj]].mu*modelgrid1[sspindx[jj]].av*$
                  (wave/5500D)^(-0.7)
                young = where(sspfits[jj].age le tbc,nyoung,comp=old,ncomp=nold)
                alam = klam*0.0
                if (nyoung ne 0) then alam[*,young] = rebin(reform(alam_bc + $
                  alam_ism,npix,1),npix,nyoung)
                if (nold ne 0) then alam[*,old] = rebin(reform(alam_ism,$
                  npix,1),npix,nold)
             endif else alam = klam*(modelgrid1[sspindx[jj]].av/rv)
             sspflux[*,*,jj] = sspflux[*,*,jj]*10.0^(-0.4*alam)
             
; do the convolution; tau=0 and no bursts is a special case
             outage = modelgrid1[sspindx[jj]].age
             if (modelgrid1[sspindx[jj]].tau eq 0.0) and $
               (modelgrid1[sspindx[jj]].nburst eq 0) then begin ; special case
                ageindx = findex(sspfits[jj].age,outage*1D9)
                outflux = interpolate(sspflux[*,*,jj],lindgen(npix),ageindx,/grid)
                outmstar = interpolate(sspfits[jj].mstar,ageindx)
                outnlyc = interpolate(sspfits[jj].nlyc,ageindx)
                outmgal = outmstar*0+1
             endif else begin
                outflux = isedfit_convolve_sfh(sspflux[*,*,jj],age=sspfits[jj].age,$
                  info=modelgrid1[sspindx[jj]],time=im_double(outage),$
                  mstar=sspfits[jj].mstar,nlyc=sspfits[jj].nlyc,$
                  cspmstar=outmstar,cspnlyc=outnlyc,nsamp=1.0,debug=debug,$
                  delayed=params.delayed,bursttype=params.bursttype)
;               test = isedfit_sfh(modelgrid1[sspindx[jj]],outage=outage,$
;                 delayed=params.delayed,bursttype=params.bursttype,/debug)
             endelse
             inf = where(finite(outflux) eq 0)
             if (inf[0] ne -1) then message, 'Bad bad bad'

; get the instantaneous and 100-Myr averaged SFRs, the birthrate
; parameter, and the total galaxy mass; normalize by OUTMGAL so that
; everything corresponds to a 1 Msun galaxy
             outsfr = isedfit_sfh(modelgrid1[sspindx[jj]],outage=outage,$
               sfr100=outsfr100,sfrage=outsfrage,b100=outb100,b_1000=outb1000,$
               mgalaxy=outmgal,delayed=params.delayed,bursttype=params.bursttype,$
               debug=debug)
             outmstar = float(outmstar/outmgal)
             outsfr = float(outsfr/outmgal)
             outsfr100 = float(outsfr100/outmgal)
             outnlyc = float(outnlyc-alog10(outmgal))
             outflux = float(outflux/outmgal)

; generate the emission-line spectrum, if desired
             if params.nebular then begin
                oiiihb = 10D^modelgrid1[sspindx[jj]].oiiihb
                nebflux = isedfit_nebular(10D^outnlyc,inst_vsigma=sspinfo.inst_vsigma,$
                  vsigma=vsigma,oiiihb=oiiihb,wave=wave,line=line,flam_line=flam_line,$
                  flam_cont=flam_cont)
; attenuate
                alam = klam[*,0]*(modelgrid1[sspindx[jj]].av/rv) ; never decrease by MU, if /CHARLOT
;               djs_plot, wave/1D4, outflux, xr=[0.1,1], xsty=3, ysty=3
;               djs_oplot, wave/1D4, outflux+nebflux, color='green'
;               djs_oplot, wave/1D4, outflux+nebflux*10.0^(-0.4*alam), color='red'
;               cc = get_kbrd(1)
                nebflux = nebflux*10.0^(-0.4*alam)
                flam_line = flam_line*10.0^(-0.4*alam)
                flam_cont = flam_cont*10.0^(-0.4*alam)

; get the EW of some emission lines; the wavelengths correspond to
; [OII], Hb+[OIII], and Ha+[NII]; should probably push this to its own
; function 
                lmin = [3727.42,4861.325,6548.043]
                lmax = [3727.42,5006.842,6730.815]
                lwave = total([[lmin],[lmax]],2)/2.0
                cflux = lwave*0.0
                for cc = 0, n_elements(lwave)-1 do begin
                   factor = sqrt(vsigma^2+sspinfo.inst_vsigma^2)/im_light(/kms)*lwave[cc]
                   cpix = where((wave lt (lmin[cc]-4*factor) and wave gt (lmin[cc]-12*factor)) or $
                     (wave gt (lmax[cc]+4*factor) and wave lt (lmax[cc]+12*factor)))
                   cflux[cc] = djs_median(outflux[cpix]+flam_cont[cpix])
                   if keyword_set(debug) then begin
                      djs_plot, wave/1D4, outflux+nebflux, xsty=3, ysty=3, $
                        xr=[lmin[cc]-30*factor,lmax[cc]+30*factor]/1D4
;                     djs_oplot, wave/1D4, outflux+flam_cont, color='red'
                      djs_oplot, wave[cpix]/1D4, outflux[cpix]+nebflux[cpix], color='green', psym=7
                      djs_oplot, [lwave[cc]]/1D4, [cflux[cc]], psym=8, color='blue', symsize=2
                      djs_oplot, !x.crange, cflux[cc]*[1,1], color='orange', line=5
                      kk = get_kbrd(1)
                   endif
                endfor
                if total(cflux le 0.0) ne 0 then message, 'This should never happen!'

                isoii = where(line.name eq '[OII]_3726' or line.name eq '[OII]_3729')
                isoiiihb = where(line.name eq '[OIII]_4959' or $
                  line.name eq '[OIII]_5007' or line.name eq 'Hbeta')
                isniiha = where(line.name eq '[NII]_6548' or $
                  line.name eq '[NII]_6584' or line.name eq 'Halpha')

                lwave = djs_mean(line[isoii].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewoii = total(line[isoii].flux*lwave)*ldust/cflux[0]

                lwave = djs_mean(line[isoiiihb].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewoiiihb = total(line[isoiiihb].flux*lwave)*ldust/cflux[1]

                lwave = djs_mean(line[isniiha].wave)
                ldust = 10.0^(-0.4*interpol(alam,wave,lwave))
                modelgrid1[sspindx[jj]].ewniiha = total(line[isniiha].flux*lwave)*ldust/cflux[2]

;               splog, modelgrid1[sspindx[jj]].age, modelgrid1[sspindx[jj]].ewoii, $
;                 modelgrid1[sspindx[jj]].ewoiiihb, modelgrid1[sspindx[jj]].ewniiha
             endif else nebflux = outflux*0.0

; ToDo: add additional interesting outputs like beta (UV slope)

; pack it in             
             modelgrid1[sspindx[jj]].mstar = outmstar
             modelgrid1[sspindx[jj]].age = outage
             modelgrid1[sspindx[jj]].sfrage = outsfrage

             modelgrid1[sspindx[jj]].sfr = alog10(outsfr>1D-15)
             modelgrid1[sspindx[jj]].sfr100 = alog10(outsfr100>1D-15)
             modelgrid1[sspindx[jj]].b100 = alog10(outb100>1D-15)
             modelgrid1[sspindx[jj]].b1000 = alog10(outb1000>1D-15)

             modelgrid1[sspindx[jj]].nlyc = outnlyc
             modelgrid1[sspindx[jj]].flux = outflux+nebflux
          endfor 
       endfor    ; close SSP loop 
       if imini eq 0 then modelgrid = temporary(modelgrid1) else $
         modelgrid = [temporary(modelgrid),temporary(modelgrid1)]
    endfor        ; close miniChunk loop   

return, modelgrid
end    

pro isedfit_montegrids, isedfit_paramfile, params=params, thissfhgrid=thissfhgrid, $
  isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, minichunksize=minichunksize, $
  priors_pdffile=priors_pdffile, debug=debug, clobber=clobber

; read the SFHGRID parameter file
    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_montegrids'
       return
    endif

    if n_elements(isedfit_dir) eq 0 then isedfit_dir = get_pwd()
    if n_elements(montegrids_dir) eq 0 then montegrids_dir = get_pwd()+'montegrids/'

; read the parameter file and then optionally call this routine
; recursively     
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)

    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
          isedfit_montegrids, params=params[ii], isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, minichunksize=minichunksize, $
            clobber=clobber, priors_pdffile=priors_pdffile, debug=debug
       endfor
       return
    endif

; number of models per output FITS table    
    modelchunksize = params.modelchunksize
    if n_elements(minichunksize) eq 0 then minichunksize = long(modelchunksize*0.1)
    
; read the SSP information structure    
    ssppath = getenv('ISEDFIT_SSP_DIR')
    if file_test(ssppath,/dir) eq 0 then begin
       splog, 'Verify that ${ISEDFIT_SSP_DIR} environment variable is defined!'
       return
    endif
    ssppath=ssppath+'/'
    sspinfofile = ssppath+'info_'+strtrim(params.spsmodels,2)+'_'+$
      strtrim(params.imf,2)+'.fits.gz'
    if file_test(sspinfofile) eq 0 then begin
       splog, 'SSP info file '+sspinfofile+' not found!'
       splog, 'Run the appropriate BUILD_*_SSP code!'
       return
    endif
    sspinfo = gz_mrdfits(sspinfofile,1,/silent)

; make directories and then check for old files
    fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,$
      montegrids_dir=montegrids_dir,priors_pdffile=priors_pdffile)

    if file_test(fp.montegrids_fullpath,/dir) eq 0 then begin
       splog, 'Making directory '+fp.montegrids_fullpath
       file_mkdir, fp.montegrids_fullpath
    endif

; ---------------------------------------------------------------------------
; first major step: build (or overwrite) the Monte Carlo grid
    montefile = strtrim(fp.montegrids_montefile,2)
    if file_test(montefile+'*') and keyword_set(clobber) eq 0 then begin
       splog, 'MonteCarlo grid file '+montefile+' exists; use /CLOBBER'
       return
    endif

; clean up old files which can conflict with this routine
    delfiles = file_search(strtrim(fp.montegrids_chunkfiles,2)+'*',count=ndel)
    if ndel ne 0 then file_delete, delfiles, /quiet

    splog, 'Building SFHGRID='+'sfhgrid'+string(params.sfhgrid,format='(I2.2)')+$
      ' REDCURVE='+strtrim(params.redcurve,2)+' IMF='+strtrim(params.imf,2)+$
      ' NMODEL='+string(params.nmodel,format='(I0)')
    montegrid = init_montegrid(params.nmodel,nmaxburst=params.nmaxburst)
       
; draw uniformly from linear TAU, or 1/TAU?
    tau = randomu(seed,params.nmodel)*(params.tau[1]-params.tau[0])+params.tau[0]
    if params.oneovertau eq 1 and params.delayed eq 1 then $
      message, 'DELAYED and ONEOVERTAU may not work well together.'
    if params.oneovertau eq 1 then tau = 1D/tau
    montegrid.tau = tau
    
; metallicity; check to make sure that the prior boundaries do not
; exceed the metallicity range available from the chosen SPSMODELS
    if (params.Zmetal[0] lt min(sspinfo.Zmetal)) then begin
       splog, 'Adjusting minimum prior metallicity!'
       params.Zmetal[0] = min(sspinfo.Zmetal)
    endif
    if (params.Zmetal[1] gt max(sspinfo.Zmetal)) then begin
       splog, 'Adjusting maximum prior metallicity!'
       params.Zmetal[1] = max(sspinfo.Zmetal)
    endif
    montegrid.Zmetal = randomu(seed,params.nmodel)*(params.Zmetal[1]-params.Zmetal[0])+params.Zmetal[0]
    
; age; unfortunately I think we have to loop to sort
    montegrid.age = randomu(seed,params.nmodel)*(params.age[1]-params.age[0])+params.age[0]
    
; reddening, if any; Gamma distribution is the default, unless FLATAV==1
    if (params.av[1] gt 0) then begin
       if params.flatav then begin
          montegrid.av = randomu(seed,params.nmodel)*(params.av[1]-params.av[0])+params.av[0] 
       endif else begin
          montegrid.av = params.av[0]*randomu(seed,params.nmodel,gamma=params.av[1])
;         montegrid.av = ((10.0^(randomn(seed,params.nmodel)*params.av[1]+alog10(params.av[0]))>0.0
;         montegrid.av = randomu(seed,params.nmodel,gamma=1.0)*params.av[1]+params.av[0]
;         montegrid.av = 10.0^(randomn(seed,params.nmodel)*params.av[1]+alog10(params.av[0]))
       endelse
       
; "mu" is the Charlot & Fall (2000) factor for evolved stellar
; populations; Gamma distribution is the default, unless FLATMU==1; do
; not let MU exceed unity
       if (strtrim(params.redcurve,2) eq 'charlot') then begin
          if params.flatmu then begin
             montegrid.mu = (randomu(seed,params.nmodel)*(params.mu[1]-params.mu[0])+params.mu[0])<1.0
          endif else begin
             montegrid.mu = (params.mu[0]*randomu(seed,params.nmodel,gamma=params.mu[1]))<1.0
;            montegrid.mu = ((10.0^(randomn(seed,params.nmodel)*params.mu[1]+alog10(params.mu[0])))<1.0)>0.0
          endelse
       endif
    endif

; add emission lines; draw [OIII]/H-beta from a uniform logarithmic
; distribution 
    if params.nebular then begin
       montegrid.oiiihb = randomu(seed,params.nmodel)*$
         (params.oiiihb[1]-params.oiiihb[0])+params.oiiihb[0] 
    endif

; now assign bursts; note that the bursts can occur outside
; (generally, before) the AGE vector; divide the time vector into
; NMAXBURST intervals of width INTERVAL_PBURST [Gyr]
    ntime = 100
;   ntime = params.nage
       
; type of burst: 0 (step function, default), 1 (gaussian), 2 (step
; function with exponential wings) 
    if (params.nmaxburst gt 0) then begin
       tburst = dblarr(params.nmaxburst,params.nmodel)-1.0
       ran = randomu(seed,params.nmaxburst,params.nmodel)
       for imod = 0L, params.nmodel-1 do begin
          for ib = 0, params.nmaxburst-1 do begin
             tmin = params.tburst[0]+ib*params.interval_pburst
             tmax = (params.tburst[0]+(ib+1)*params.interval_pburst)<params.tburst[1]
;            tmin = params.minage+ib*params.interval_pburst
;            tmax = (params.minage+(ib+1)*params.interval_pburst)<params.maxage

             if (tmax le tmin) then message, 'This violates causality!'
             time1 = randomu(seed,ntime)*(tmax-tmin)+tmin
             time1 = time1[sort(time1)]
                
; compute the cumulative probability, dealing with edge effects correctly
             dtimeleft = (time1-shift(time1,1))/2.0
             dtimeleft[0] = time1[0]-tmin
             dtimeright = (shift(time1,-1)-time1)/2.0
             dtimeright[ntime-1] = tmax-time1[ntime-1]
             prob = params.pburst*total([[dtimeleft],[dtimeright]],2)/$
               params.interval_pburst
             if (prob[0] le 0) then message, 'This should not happen!'
             
             this = where(total(prob,/cum) gt ran[ib,imod])
             if (this[0] ne -1) then tburst[ib,imod] = time1[this[0]]
;            splog, tmin, tmax, ran[ib,imod], tburst[ib,imod] & if ib eq 0 then print
          endfor
       endfor 
       montegrid.nburst = total(tburst gt -1.0,1) ; total number of bursts
    endif
    
; assign burst strengths and durations       
    hasburst = where(montegrid.nburst gt 0,nhasburst)
    nallburst = long(total(montegrid.nburst))
    if (nhasburst ne 0L) then begin
       if params.flatdtburst then begin
          dtburst = randomu(seed,nallburst)*(params.dtburst[1]-$
            params.dtburst[0])+params.dtburst[0]
       endif else begin
          dtburst = 10.0^(randomu(seed,nallburst)*(alog10(params.dtburst[1])-$
            alog10(params.dtburst[0]))+alog10(params.dtburst[0]))
       endelse
       
       if params.flatfburst then begin
          fburst = randomu(seed,nallburst)*(params.fburst[1]-$
            params.fburst[0])+params.fburst[0]
       endif else begin
          fburst = 10.0^(randomu(seed,nallburst)*(alog10(params.fburst[1])-$
            alog10(params.fburst[0]))+alog10(params.fburst[0]))
       endelse
       
       count = 0L
       for ii = 0L, nhasburst-1 do begin ; sort
          if (montegrid[hasburst[ii]].nburst gt 0) then begin
             nb = montegrid[hasburst[ii]].nburst
             good = where(tburst[*,hasburst[ii]] gt -1.0)
             montegrid[hasburst[ii]].tburst = tburst[good,hasburst[ii]]
             montegrid[hasburst[ii]].fburst = fburst[count:count+nb-1]
             montegrid[hasburst[ii]].dtburst = dtburst[count:count+nb-1]
             count = count + nb
          endif
       endfor
       
; allow a FRACTRUNC fraction of the models with bursts to have
; truncated bursts 
       if (params.fractrunc gt 0D) then begin
          ntrunc = long(params.fractrunc*nhasburst)
          if params.fractrunc eq 1D then trunc = lindgen(ntrunc) else $
            trunc = random_indices(nhasburst,ntrunc)
       endif else ntrunc = 0L
       if (ntrunc gt 0L) then begin
          if (params.bursttype ne 1) then message, 'Should use truncated *Gaussian* bursts.'
;         montegrid[hasburst[trunc]].trunctau = params.trunctau ; [Gyr]
          montegrid[hasburst[trunc]].trunctau = $
            randomu(seed,ntrunc)*(params.trunctau[1]-$
            params.trunctau[0])+params.trunctau[0] ; [Gyr]
       endif
    endif  
    
; ---------------------------------------------------------------------------
; now build the actual composite stellar populations; do it in chunks
; to avoid memory issues and pass along all the parameters

; loop on each chunk of models
    t0 = systime(1)
    nchunk = params.nmodelchunk
    chunkfiles = strtrim(fp.montegrids_chunkfiles,2)
    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*modelchunksize
       i2 = ((ichunk*modelchunksize+modelchunksize)<params.nmodel)-1L
       nthese = i2-i1+1L
       these = lindgen(nthese)+i1

       modelgrid = build_modelgrid(montegrid[these],params=params,debug=debug,$
         sspinfo=sspinfo,ssppath=ssppath+strtrim(params.spsmodels,2)+'/',$
         minichunksize=minichunksize,ichunk=ichunk,nchunk=nchunk)
       modelgrid.chunkindx = ichunk

       im_mwrfits, modelgrid, chunkfiles[ichunk], /clobber

; build an updated montegrid structure which also contains the derived
; parameters we want to look at
       newmontegrid1 = struct_trimtags(modelgrid,except=['wave','flux'])
       if ichunk eq 0 then newmontegrid = temporary(newmontegrid1) else $
         newmontegrid = [temporary(newmontegrid),temporary(newmontegrid1)]
    endfor
    splog, 'Total time (min) = ', (systime(1)-t0)/60.0

; write out MONTEFILE
    im_mwrfits, newmontegrid, montefile, /clobber

; finally, build a QAplot; allow the user to overwrite PDFFILE 
    if file_test(fp.isedfit_dir+fp.qaplot_priors_pdffile) and $
      keyword_set(clobber) eq 0 then begin
       splog, 'Output file '+fp.qaplot_priors_pdffile+' exists; use /CLOBBER'
       return
    endif

    qafile = fp.isedfit_dir+fp.qaplot_priors_psfile
    montegrids_qaplot, newmontegrid, params=params, qafile=qafile

return
end
