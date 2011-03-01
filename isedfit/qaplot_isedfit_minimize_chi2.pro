pro qaplot_isedfit_minimize_chi2, xvalue, pdf, quant=quant, $
  best=best, median=median, chi2=chi2, type=type
; build a simple QAplot for inspection (called by
; ISEDFIT_MINIMIZE_CHI2) 

    im_plotconfig, 6, pos
    xrange = minmax(xvalue)
    im_plothist, xvalue, xx, yy, weight=pdf, bin=0.04, /noplot
    djs_plot, xx, yy, position=pos[*,0], xsty=3, ysty=3, $
      xrange=xrange, color=djs_icolor('red'), $
      xtickname=replicate(' ',10), psym=10
    djs_oplot, quant[0]*[1,1], !y.crange, line=2, thick=2
    djs_oplot, quant[2]*[1,1], !y.crange, color='cyan'
    djs_oplot, quant[5]*[1,1], !y.crange, color='cyan'

    label = [type+'_{best} = '+strtrim(string(best,format='(F12.3)'),2),$
     type+'_{50} = '+strtrim(string(median,format='(F12.3)'),2),$
     '\chi^{2}_{min} = '+strtrim(string(chi2,format='(F12.1)'),2)]
    im_legend, label, /left, /top, box=0, margin=0, charsize=1.6

    srt = sort(xvalue)
    djs_plot, xvalue[srt], total(pdf[srt],/cum,/double), /noerase, $
      position=pos[*,1], xsty=3, ysty=3, xrange=xrange, $
      color='red', xtitle=type, yrange=[0,1]
    djs_oplot, quant[0]*[1,1], !y.crange, line=2, thick=2
    djs_oplot, quant[2]*[1,1], !y.crange, color='cyan'
    djs_oplot, quant[5]*[1,1], !y.crange, color='cyan'

return
end    
    
