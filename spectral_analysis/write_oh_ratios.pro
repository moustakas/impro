pro write_oh_ratios
; jm04jul11uofa
; write out Gaussian look-up tables for the emission-line ratios that
; are most sensitive to electron-temperature abundance

    hii = read_hii_regions()
    good = where((hii.ZT_log12oh gt -900.0) and (hii.oiii_5007_oii gt -900.0) and $
      (hii.oii_nii_6584 gt -900.0) and (hii.nii_6584_h_alpha gt -900.0),ngood)
    hii = hii[good]

    keeptags = ['nii_6584_h_alpha','nii_6584_h_alpha_err',$
      'nii_6584_oii','nii_6584_oii_err','oiii_5007_oii','oiii_5007_oii_err']
    ratios = ['nii_ha','nii_ha_err','nii_oii','nii_oii_err','oiii_oii','oiii_oii_err']
    labels = ['[N II]/Ha','[N II]/[O II]','[O III]/[O II]']
    newtags = [ratios]
    info = im_struct_trimtags(hii,select=keeptags,newtags=newtags)

    oh = hii.zt_log12oh
    oh_err = hii.zt_log12oh_err

    nratio = n_elements(labels)

    result = {$
      label:     '', $
      coeff:     fltarr(2), $
      coeff_err: fltarr(2)}
    result = replicate(result,nratio)

    maxiter = 5L
    nsigma = 3L
    ndegree = 1L
    
    for i = 0L, nratio-1L do begin

       x = info.(2*i)
       xerr = info.(2*i+1L)

       srt = sort(x)
       x = x[srt]
       xerr = xerr[srt]

       y = oh[srt]
       yerr = oh_err[srt]

       good = bytarr(n_elements(x))+1B
       w = lindgen(n_elements(x))
       
       for iter = 0L, maxiter-1L do begin
          
          coeff = poly_fit(x[w],y[w],ndegree,measure_errors=yerr[w],$
            sigma=sigma,yfit=yfit,/double)
          res = y[w] - yfit
          sig = stddev(res)
          good[w] = good[w]*(abs(res) le nsigma*sig)
          w = where(good)

       endfor

       rejected = where(good eq 0B,nrejected,comp=kept)
       
       xaxis = findgen((max(x)-min(x))/0.01+1)*0.01+min(x)
       yfit = poly(xaxis,coeff)
       
       ploterror, x[kept], y[kept], xerr[kept], yerr[kept], $
         ps=4, xsty=3, ysty=3, charsize=2.0, charthick=2.0
       oploterror, x[rejected], y[rejected], xerr[rejected], yerr[rejected], $
         ps=4, errcolor=djs_icolor('yellow'), color=djs_icolor('yellow')
       djs_oplot, xaxis, yfit, line=0, thick=2.0, color='red'
       cc = get_kbrd(1)
       
    endfor    

stop    
    
; ---------------------------------------------------------------------------
; Danielle's Gaussian sampling    
; ---------------------------------------------------------------------------
    
    bin = 0.25
    minpoints = 30L
    nsamp = 10.0

    nratio = n_elements(ratios)
    for i = 0L, nratio-1L do begin

       y = oh
       x = info.(2*i)
       srt = sort(x)
       x = x[srt]
       y = y[srt]

       xmin = min(x)
       xmax = max(x)

       xindx = lindgen(ceil((xmax-xmin)/bin+1))*bin+xmin
       nbins = n_elements(xindx)

       for j = 0L, nbins-2L do begin

          ninrange = 0L
          count = 0.0
          okay = 0L

          while (ninrange lt minpoints) and (okay eq 0L) do begin

             xlo = xindx[j]
             xhi = xindx[j+1L]+count*bin/4.0
             
             inrange = where((x gt xlo) and (x lt xhi),ninrange)
             if (xhi gt xmax) then okay = 1L

             count = count + 1.0
             
          endwhile

; reject outliers

          djs_iterstat, y[inrange], sigrej=6.0, mask=mask
          keep = where(mask eq 1B,nkeep)
          
          plothist, y[inrange[keep]], xhist, yhist, /noplot, bin=0.05
          
          ygauss = mpfitpeak(xhist,yhist,a,/gauss,/positive,nterms=3)
          xgauss = (findgen(nsamp*((max(xhist)-min(xhist))/0.05))*0.05/nsamp + min(xhist))
          ygauss = mpfitpeak_gauss(xgauss,a)

          if (n_elements(xmean) eq 0L) then xmean = djs_mean(x[inrange]) else $
            xmean = [xmean,djs_mean(x[inrange])]
          if (n_elements(gmean) eq 0L) then gmean = a[1] else gmean = [gmean,a[1]]
          if (n_elements(gsigma) eq 0L) then gsigma = a[2] else gsigma = [gsigma,a[2]]
          
          djs_plot, xhist, yhist, ps=10, xsty=3, ysty=3, xthick=2.0, $
            ythick=2.0, thick=2.0
          djs_oplot, xgauss, ygauss, line=0, thick=2.0, color='red', ps=10
;         cc = get_kbrd(1)
          
       endfor

       gbig = interpol(gmean,xmean,x)
       gsigmabig = interpol(gsigma,xmean,x)
       
       djs_plot, x, y, ps=4, xsty=3, ysty=3
;      djs_oplot, x, gbig, line=0, thick=2.0, color='red'
       oploterror, x, gbig, gsigmabig, line=0, thick=2.0, $
         color=djs_icolor('red'), errcolor=djs_icolor('red')
       
stop
       
    endfor
    
    stop

return
end
   
