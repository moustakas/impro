function manygauss, xindx, pp, nline=nline, nback=nback, loglam=loglam, $
  sigmares=sigmares, background=background

    yval = 0.0D

    for iline=0, nline-1 do $
      yval = yval + onegauss(loglam[xindx], pp[iline*3:iline*3+2], sigmares=sigmares[iline])
    for iback=0, nback-1 do $
      yval = yval + background[xindx,iback] * pp[nline*3+iback]

;   plot, xindx, yval
;   cc = get_kbrd(1)
    
return, yval
end

