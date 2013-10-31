;+
; NAME:
;   ISEDFIT_QAPLOT_MODELS
;
; PURPOSE:
;   Generate color-color and redshift-color QAplots to help assess the
;   choice of iSEDfit model parameters/priors. 
;
; INPUTS:
;   isedfit_paramfile - iSEDfit parameter file
;   maggies - input galaxy photometry [NFILT,NGAL]
;   ivarmaggies - inverse variance array for MAGGIES [NFILT,NGAL]  
;   z - input galaxy redshifts [NGAL] 
;
; OPTIONAL INPUTS:
;   params - data structure with the same information contained in
;     ISEDFIT_PARAMFILE (over-rides ISEDFIT_PARAMFILE)
;   thissfhgrid - if ISEDFIT_PARAMFILE contains multiple grids then
;     build this SFHgrid (may be a vector)
;   isedfit_dir - full directory path where the QAplots should be
;     written; must match the directory passed to ISEDFIT_MODELS
;     (default PWD=present working directory) 
;   thesefilters - build the QAplots from this subset of the available
;     filters; default is to use all the filters, which is *not*
;     recommended if you have more than 4-5 filters (see COMMENTS)!
;     the (string) filter names can be missing the '.par' suffix
;   outprefix - optional output prefix string (see ISEDFIT) 
;
;   colorcolor_pdffile - overwrite the default name of the output
;     color-color PDF file (not typically necessary since the file
;     name matches the ISEDFIT output files); must end with a '.PDF'
;     suffix; if COLORCOLOR_PDFFILE='' then don't make this plot
;   zcolor_pdffile - same as above but for the redshift-color PDF
;     file; if ZCOLOR_PDFFILE='' then don't make this plot
;
; KEYWORD PARAMETERS:
;   clobber - overwrite existing files of the same name (the default
;     is to check for existing files and if they exist to exit
;     gracefully)  
;
; OUTPUTS:
;   This routine generates two sets of QAplots--color-color diagrams
;   and color-redshift diagrams---which show how well (or poorly!) the
;   models overlap the data.  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Left to its own devices, this routine will generate *all*
;   combinations of color-redshift and color-color diagrams.
;   Therefore if you have more than a handful of filters you should
;   use the THESEFILTERS optional input parameter.  For example,
;   consider your full list of filters consists of:
;       filterlist = [$
;         'galex_FUV.par',$
;         'galex_NUV.par',$
;         'sdss_u0.par',$
;         'sdss_g0.par',$
;         'sdss_r0.par',$
;         'sdss_i0.par',$
;         'sdss_z0.par',$
;         'wise_w1.par',$
;         'wise_w2.par']
;
;   You might call this routine with the following (intelligently
;   chosen) subset of filters:
;     IDL> thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
;     IDL> isedfit_qaplot_models, isedfit_paramfile, maggies, $
;     IDL>  ivarmaggies, z, thesefilters=thesefilters
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Aug 19, Siena 
;
; Copyright (C) 2013, John Moustakas
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
  thesefilters=thesefilters, outprefix=outprefix, colorcolor_pdffile=colorcolor_pdffile, $
  zcolor_pdffile=zcolor_pdffile, clobber=clobber

    if n_elements(isedfit_paramfile) eq 0 and n_elements(params) eq 0 then begin
       doc_library, 'isedfit_qaplot_models'
       return
    endif

; read the parameter file; parse to get the relevant path and
; filenames
    if (n_elements(params) eq 0) then params = $
      read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
    if (n_elements(isedfit_dir) eq 0) then isedfit_dir = get_pwd()

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
          isedfit_qaplot_models, isedfit_paramfile1, maggies1, ivarmaggies1, z, $
            params=params[ii], isedfit_dir=isedfit_dir, thesefilters=thesefilters, $
            outprefix=outprefix, colorcolor_pdffile=colorcolor_pdffile, $
            zcolor_pdffile=zcolor_pdffile, clobber=clobber
       endfor
       return
    endif

; allow the user to overwrite PDFFILE
    fp = isedfit_filepaths(params,outprefix=outprefix,isedfit_dir=isedfit_dir,$
      colorcolor_pdffile=colorcolor_pdffile,zcolor_pdffile=zcolor_pdffile)

    if fp.qaplot_zcolor_pdffile ne '' then begin
       if file_test(fp.isedfit_dir+fp.qaplot_zcolor_pdffile) and $
         keyword_set(clobber) eq 0 then begin
          splog, 'Output file '+fp.qaplot_zcolor_pdffile+' exists; use /CLOBBER'
          return
       endif
    endif
    
    if fp.qaplot_colorcolor_pdffile ne '' then begin
       if file_test(fp.isedfit_dir+fp.qaplot_colorcolor_pdffile) and $
         keyword_set(clobber) eq 0 then begin
          splog, 'Output file '+fp.qaplot_colorcolor_pdffile+' exists; use /CLOBBER'
          return
       endif
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
    if params.nzz eq 1 then zindx = 0 else $
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

    if params.nzz eq 1 then modelmaggies = models.modelmaggies[filtmatch,*] else $
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
    npix = ((long(nmodel*1E-3)+(odd(long(nmodel*1E-3)) eq 0))>21L)<51L
    mylevels = [0.5,0.75,0.99]
    mycann = string(mylevels,format='(F4.2)')
    mylevels = [0.5,0.75,0.99]
    
    nperpage = 4
    
; make the redshift-color plots
    npage = ceil(ncombos/float(nperpage))

    if fp.qaplot_zcolor_pdffile ne '' then begin
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
                if ngood gt 500 then begin
                   im_hogg_scatterplot, z[good], ygal, /overplot, $
                     /outlier, outpsym=symcat(6,thick=1), outsymsize=0.3, /nogrey, $
                     outcolor=im_color('dodger blue',100), contour_color=im_color('navy',101), $
                     xrange=xrange, yrange=yrange, position=pos[*,ii], xsty=7, ysty=5, $
                     levels=mylevels, /internal, c_annotation=mycann
                endif else begin
                   djs_oplot, z[good], ygal, psym=symcat(6,thick=3), symsize=0.5, $
                     color=im_color('dodger blue',100)
                endelse
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
;               oplot, bigredshift, ymodel, psym=3, color=im_color('red')
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
    endif

; make all the color-color plot combinations (ignore redshift for now) 
    npage = ceil(n_elements(ccinfo)/float(nperpage))

    if fp.qaplot_colorcolor_pdffile ne '' then begin
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
                
                contour_color = [im_color('red',252),im_color('forest green',253),$
                  im_color('grey40',254)]
                im_hogg_scatterplot, xmodel, ymodel, position=pos[*,ii], xsty=1, ysty=1, $
                  xrange=xrange, yrange=yrange, xtitle=ccinfo[count].xtitle, levels=mylevels, $
                  ytitle=ccinfo[count].ytitle, noerase=ii gt 0, outlier=0, /nogrey, $
                  contour_color=[252,253,254], cthick=6, $
                  xnpix=npix, ynpix=npix, /internal, c_annotation=mycann
                if ngood gt 500 then begin
                   im_hogg_scatterplot, xgal, ygal, /overplot, $
                     /outlier, outpsym=symcat(6,thick=1), outsymsize=0.3, /nogrey, $
                     outcolor=im_color('dodger blue',100), contour_color=im_color('navy',101), $
                     levels=mylevels, xrange=xrange, yrange=yrange, position=pos[*,ii], $
                     xnpix=npix, ynpix=npix, /internal, c_annotation=mycann
                endif else begin
                   djs_oplot, xgal, ygal, psym=symcat(6,thick=3), symsize=0.5, $
                     color=im_color('dodger blue',100)
                endelse
                count++
             endif
          endfor
          xyouts, pos[0,0], pos[3,0]+0.14, 'Model Photometry:', align=0.0, $
            /normal, charsize=1.3
          im_legend, ['50 Percentile','75 Percentile','98 Percentile'], $
            box=1, /left, /top, line=[0,5,3], pspacing=1.9, $ ; psym=-[16,15,14], 
            position=[pos[0,0],pos[3,0]+0.12], /norm, charsize=1.3, $
            color=['red','forest green','grey40'], thick=6 ;, symsize=[1.2,1.2,1.4]
          
          im_legend, 'Galaxy Photometry', box=1, /left, /top, $
            line=0, pspacing=1.9, position=[pos[0,1],pos[3,1]+0.08], /norm, $
            charsize=1.3, color='navy', thick=6
          plots, pos[0,1]+0.055, pos[3,1]+0.062, psym=-symcat(6,thick=6), $
            symsize=1.5, color=im_color('dodger blue'), /norm
       endfor 
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

return
end
