;+
; NAME:
;   ISEDFIT_QAPLOT_MODELS
;
; PURPOSE:
;   Generate color-redshift quality-assurance (QA) plots from the
;   iSEDfit output.
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

pro isedfit_qaplot_models, isedfit_paramfile, maggies1, ivarmaggies1, z, $
  params=params, thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, $
  montegrids_dir=montegrids_dir, thesefilters=thesefilters, $
  pdffile=pdffile, clobber=clobber

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_qaplot_models'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()
    if (n_elements(montegrids_dir) eq 0) then montegrids_dir = get_pwd()+'montegrids/'

; error checking on the input photometry    
    ndim = size(maggies1,/n_dim)
    dims = size(maggies1,/dim)
    if (ndim eq 1) then ngal = 1 else ngal = dims[1]  ; number of galaxies
    nfilt = dims[0] ; number of filters

    nmaggies = n_elements(maggies1)
    nivarmaggies = n_elements(ivarmaggies1)
    nz = n_elements(z)
    
    if (nmaggies eq 0L) or (nivarmaggies eq 0L) or $
      (nz eq 0L) then begin
       doc_library, 'isedfit'
       return
    endif

    ndim = size(maggies1,/n_dimension)
    if (ndim ne 2L) then begin ; reform into a 2D array
       maggies1 = reform(maggies1,n_elements(maggies1),1)
       ivarmaggies1 = reform(ivarmaggies1,n_elements(maggies1),1)
    endif

    if (n_elements(maggies1) ne n_elements(ivarmaggies1)) then begin
       splog, 'Dimensions of MAGGIES and IVARMAGGIES do not match'
       return
    endif
    if (nz ne ngal) then begin
       splog, 'Dimensions of MAGGIES and Z do not match'
       return
    endif
    if (total(finite(maggies1) eq 0) ne 0.0) or $
      (total(finite(ivarmaggies1) eq 0) ne 0.0) then begin
       splog, 'MAGGIES and IVARMAGGIES cannot have infinite values!'
       return
    endif
    if (total(z le 0.0) ne 0.0) then begin
       splog, 'Z should all be positive'
       return
    endif

; treat each SFHgrid separately
    ngrid = n_elements(params)
    if ngrid gt 1 then begin
       for ii = 0, ngrid-1 do begin
stop 
          
          isedfit_qaplot_models, params=params[ii], isedfit_dir=isedfit_dir, $
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
    if file_test(fp.isedfit_dir+fp.qaplot_zcolor_pdffile) and $
      keyword_set(clobber) eq 0 then begin
       splog, 'Output file '+fp.qaplot_zcolor_pdffile+' exists; use /CLOBBER'
       return
    endif
    
; get the (desired) filters and the redshift grid
    filterlist = strtrim(params.filterlist,2)
;   filterlist = ['U','B','V','R','I','Y'] ; test

    if n_elements(thesefilters) eq 0 then thesefilters = filterlist
    if size(thesefilters,/type) ne 7 then begin
       splog, 'THESEFILTERS must be type string!'
       return
    endif
    filtmatch = intarr(n_elements(thesefilters))-1
    for ii = 0, n_elements(thesefilters)-1 do filtmatch[ii] = $
      where(strmatch(filterlist,'*'+thesefilters[ii]+'*'))
    check = where(filtmatch eq -1)
    if check[0] ne -1 then begin
       splog, 'No match for these THESEFILTERS: '+strjoin(thesefilters[check],', ')
       return
    endif

    filterlist = filterlist[filtmatch]
    maggies = maggies1[filtmatch,*]
    ivarmaggies = ivarmaggies1[filtmatch,*]
    
    short = '['+repstr(filterlist,'.par','')+']'
    nfilt = n_elements(filterlist)
    ncombos = nfilt*(nfilt-1)/2 ; number of combinations

    splog, 'Building QAplots using the following filters:'
    niceprint, replicate('  ',nfilt), filterlist
    
    redshift = params.redshift
    nredshift = n_elements(redshift)
    zindx = findex(params.redshift,redshift)
;   zindx = findex(redshift,z) ; test code

; restore all the model photometry
    chunkfile = strtrim(fp.models_chunkfiles[0],2) ; check the first file
    if file_test(chunkfile+'.gz') eq 0 then begin
       splog, 'First ChunkFile '+chunkfile+' not found!'
       return
    endif

    splog, 'Restoring the model photometry'
    nchunk = n_elements(fp.models_chunkfiles)
    for ichunk = 0, nchunk-1 do begin
       splog, 'Reading '+fp.models_chunkfiles[ichunk]+'.gz'
       models1 = mrdfits(fp.models_chunkfiles[ichunk]+'.gz',1,/silent)
       if ichunk eq 0L then models = temporary(models1) else $
         models = [temporary(models),temporary(models1)]
    endfor
    nmodel = n_elements(models)
    if nmodel ne params.nmodel then message, 'Mismatched number of models!'

    modelmaggies = interpolate(transpose(models.modelmaggies[filtmatch,*],[0,2,1]),zindx)
    bigredshift = rebin(reform(redshift,1,nredshift),nmodel,nredshift)

; compute all the possible combinations of model and galaxy magnitudes 
    colortitle = strarr(ncombos)
    galaxymag = fltarr(ncombos,ngal)-1.0
    modelmag = fltarr(ncombos,nmodel,nredshift)

    splog, 'Computing AB magnitudes for the model and galaxy samples.'
    count = 0
    for jj = 0, nfilt-2 do for kk = jj+1, nfilt-1 do begin
       colortitle[count] = short[jj]+'-'+short[kk]
; check for no photometry
       good = where(maggies[jj,*] gt 0.0 and maggies[kk,*] gt 0.0,ngood)
       if ngood ne 0L then galaxymag[count,good] = $
         -2.5*alog10(maggies[jj,good]/maggies[kk,good])
       modelmag[count,*,*] = -2.5*alog10(modelmaggies[jj,*,*]/modelmaggies[kk,*,*])
       count++
    endfor

; this little algorithm builds all the color-color combinations we
; will need for the plot below
;   splog, 'Building color-color combinations.'
    ccinfo1 = {xgal: fltarr(ngal), ygal: fltarr(ngal), xmodel: fltarr(nmodel,nredshift), $
      ymodel: fltarr(nmodel,nredshift), xtitle: '', ytitle: ''}
    for bluecol = 0, nfilt-2 do begin ; blue color
       if bluecol eq 0 then bluestart = 0 else bluestart = ib
       nblue = nfilt-1-bluecol  ; lose one row per column
       bcount = 1
       for ib = bluestart, bluestart+nblue-1 do begin
; ignore the bluest and reddest filter combination       
          if bluecol eq 0 and ib eq nfilt-1 then continue
          nred = nblue-bcount
          if ib eq bluestart then redstart = nblue+ib
          rcount = 0
          for ir = redstart, redstart+nred-1 do begin
             if n_elements(ccinfo) eq 0L then ccinfo = ccinfo1 else $
               ccinfo = [temporary(ccinfo),ccinfo1]
             indx = n_elements(ccinfo)-1
             ccinfo[indx].xgal = galaxymag[ir,*]
             ccinfo[indx].ygal = galaxymag[ib,*]
             ccinfo[indx].xmodel = modelmag[ir,*,*]
             ccinfo[indx].ymodel = modelmag[ib,*,*]
             ccinfo[indx].xtitle = colortitle[ir]
             ccinfo[indx].ytitle = colortitle[ib]
;            print, nblue, nred, ib, ir, ' ', colortitle[ib], ' ', colortitle[ir]
             rcount++
             nred = nred-rcount
          endfor
          bcount++
          redstart = ir
          bluestart += nblue
       endfor 
    endfor

; some plotting preferences
    npix = ((long(nmodel*1E-3)+(odd(long(nmodel*1E-3)) eq 0))>21L)<101L
    mylevels = [0.5,0.75,0.98]
    mycann = string(mylevels,format='(F4.2)')
    mylevels = [0.5,0.75,0.98]
    
    nperpage = 4
    
; make the redshift-color plots
    npage = ceil(ncombos/float(nperpage))

    psfile = fp.isedfit_dir+fp.qaplot_zcolor_psfile
    im_plotconfig, 5, pos, psfile=psfile, xspace=1.0, $
      yspace=0.8, xmargin=[1.0,0.6], width=3.05*[1,1], $
      height=2.6*[1,1], charsize=1.4, ymargin=[2.0,1.1]

    xrange = minmax(redshift)
    quant = fltarr(nredshift,6) ; quantiles

    count = 0
    for pp = 0, npage-1 do begin
       for ii = 0, nperpage-1 do begin
          if count le ncombos-1 then begin
             ymodel = reform(modelmag[count,*,*])
             for jj = 0, nredshift-1 do quant[jj,*] = weighted_quantile($
               ymodel[*,jj],quant=[0.25,0.75,0.05,0.95,0.0,1.0])
             
             good = where(galaxymag[count,*] gt -1.0,ngood)
             if ngood ne 0L then begin
                ygal = reform(galaxymag[count,good])
                yrange = [min(ymodel)<weighted_quantile(ygal,quant=0.01),$
                  max(ymodel)>weighted_quantile(ygal,quant=0.99)]
             endif else begin
                yrange = minmax(ymodel)
             endelse
             yrange += (yrange[1]-yrange[0])*0.05*[-1,1]
             
             plot, [0], [0], /nodata, position=pos[*,ii], xsty=3, ysty=1, $
               xrange=xrange, yrange=yrange, xtitle='Redshift', $
               ytitle=colortitle[count], noerase=ii gt 0, xticks=3
             im_hogg_scatterplot, z[good], ygal, /overplot, $
               /outlier, outpsym=symcat(6,thick=1), outsymsize=0.3, /nogrey, $
               outcolor=im_color('dodger blue',100), contour_color=im_color('navy',101), $
               xrange=xrange, yrange=yrange, position=pos[*,ii], xsty=7, ysty=5, $
               levels=mylevels, /internal, c_annotation=mycann
;            if ngood ne 0L then oplot, z[good], ygal, psym=symcat(6,thick=1), $
;              color=im_color('blue'), symsize=0.2
             oplot, redshift, quant[*,0], line=0, color=im_color('red'), $
               thick=4, psym=-symcat(16), symsize=0.7
             oplot, redshift, quant[*,1], line=0, color=im_color('red'), $
               thick=4, psym=-symcat(16), symsize=0.7
             oplot, redshift, quant[*,2], line=5, color=im_color('forest green'), $
               thick=4, psym=-symcat(15), symsize=0.7
             oplot, redshift, quant[*,3], line=5, color=im_color('forest green'), $
               thick=4, psym=-symcat(15), symsize=0.7
             oplot, redshift, quant[*,4], line=3, color=im_color('grey40'), $
               thick=4, psym=-symcat(14), symsize=0.7
             oplot, redshift, quant[*,5], line=3, color=im_color('grey40'), $
               thick=4, psym=-symcat(14), symsize=0.7
;            oplot, bigredshift, ymodel, psym=3, color=im_color('red')
             count++
          endif
       endfor 
       xyouts, pos[0,0], pos[3,0]+0.14, 'Model Photometry:', align=0.0, $
         /normal, charsize=1.3
       im_legend, ['25%, 75% Quantile',' 5%, 95% Quantile',' 0%,100% Quantile'], $
         box=1, /left, /top, line=[0,5,3], psym=-[16,15,14], pspacing=1.9, $
         position=[pos[0,0],pos[3,0]+0.12], /norm, charsize=1.3, $
         color=['red','forest green','grey40'], thick=6, symsize=[1.2,1.2,1.4]

       im_legend, 'Galaxy Photometry', box=1, /left, /top, $
         line=0, pspacing=1.9, position=[pos[0,1],pos[3,1]+0.08], /norm, $
         charsize=1.3, color='navy', thick=6
       plots, pos[0,1]+0.055, pos[3,1]+0.062, psym=-symcat(6,thick=6), $
         symsize=1.5, color=im_color('dodger blue'), /norm
    endfor 

    im_plotconfig, psfile=psfile, /psclose, /pdf

; make all the color-color plot combinations (ignore redshift for now) 
    npage = ceil(n_elements(ccinfo)/float(nperpage))

    psfile = fp.isedfit_dir+fp.qaplot_colorcolor_psfile
    im_plotconfig, 5, pos, psfile=psfile, xspace=1.0, $
      yspace=0.8, xmargin=[1.0,0.6], width=3.05*[1,1], $
      height=2.6*[1,1], charsize=1.4, ymargin=[2.0,1.1]

    count = 0
    for pp = 0, npage-1 do begin
       for ii = 0, nperpage-1 do begin
          if count lt n_elements(ccinfo) then begin
             good = where(ccinfo[count].xgal gt -1.0 and ccinfo[count].ygal gt -1.0,ngood)
             xgal = reform(ccinfo[count].xgal[good])
             ygal = reform(ccinfo[count].ygal[good])
             xmodel = reform(ccinfo[count].xmodel)
             ymodel = reform(ccinfo[count].ymodel)

             yrange = [min(ymodel)<weighted_quantile(ygal,quant=0.01),$
               max(ymodel)>weighted_quantile(ygal,quant=0.99)]
             xrange = [min(xmodel)<weighted_quantile(xgal,quant=0.01),$
               max(xmodel)>weighted_quantile(xgal,quant=0.99)]
             
             contour_color = [im_color('red',252),im_color('forest green',253),im_color('grey40',254)]
             im_hogg_scatterplot, xmodel, ymodel, position=pos[*,ii], xsty=1, ysty=1, $
               xrange=xrange, yrange=yrange, xtitle=ccinfo[count].xtitle, levels=mylevels, $
               ytitle=ccinfo[count].ytitle, noerase=ii gt 0, outlier=0, /nogrey, $
               contour_color=[252,253,254], cthick=6, $
               xnpix=npix, ynpix=npix, /internal, c_annotation=mycann
             im_hogg_scatterplot, xgal, ygal, /overplot, $
               /outlier, outpsym=symcat(6,thick=1), outsymsize=0.3, /nogrey, $
               outcolor=im_color('dodger blue',100), contour_color=im_color('navy',101), $
               levels=mylevels, xrange=xrange, yrange=yrange, position=pos[*,ii], $
               xnpix=npix, ynpix=npix, /internal, c_annotation=mycann
;            djs_oplot, xgal, ygal, psym=3, $ ; psym=symcat(8,thick=1), symsize=0.5, $
             count++
          endif
       endfor
       xyouts, pos[0,0], pos[3,0]+0.14, 'Model Photometry:', align=0.0, $
         /normal, charsize=1.3
       im_legend, ['50 Percentile','75 Percentile','98 Percentile'], $
         box=1, /left, /top, line=[0,5,3], pspacing=1.9, $ ; psym=-[16,15,14], 
         position=[pos[0,0],pos[3,0]+0.12], /norm, charsize=1.3, $
         color=['red','forest green','grey40'], thick=6;, symsize=[1.2,1.2,1.4]

       im_legend, 'Galaxy Photometry', box=1, /left, /top, $
         line=0, pspacing=1.9, position=[pos[0,1],pos[3,1]+0.08], /norm, $
         charsize=1.3, color='navy', thick=6
       plots, pos[0,1]+0.055, pos[3,1]+0.062, psym=-symcat(6,thick=6), $
         symsize=1.5, color=im_color('dodger blue'), /norm
    endfor 
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
