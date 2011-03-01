function oh2r23, logoh, logq=logq1, o32=o32, debug=debug
; jm06apr18uofa - given an oxygen abundance, return the R23 value
;                 (so far only supports the KK04 calibration)

    nlogoh = n_elements(logoh)
    if (nlogoh eq 0L) then return, -1L

    nlogq = n_elements(logq1)
    case nlogq of
       0L: logq = replicate(alog10(4D7),nlogoh) ; note!
       1L: logq = replicate(logq1,nlogoh)
       else: logq = logq1
    endcase
    nlogq = n_elements(logq)

    if (nlogq ne nlogoh) then begin
       splog, 'DIMENSIONS of LOGQ and LOGOH must match.'
       return, -1L
    endif
    
    if (nlogoh gt 1L) then begin
       r23 = fltarr(nlogoh) & o32 = fltarr(nlogoh)
       for i = 0L, nlogoh-1L do begin
          r23[i] = oh2r23(logoh[i],logq=logq[i],o32=o321,debug=debug)
          o32[i] = o321
       endfor
       return, r23
    endif
    
    if (logq lt alog10(5D6)) or (logq gt alog10(3D8)) then begin
       splog, 'LOGQ outside model boundaries!'
       return, -1L
    endif

    maxr23 = 1.1 & minr23 = -1.0 & dr23 = 0.01
    r23axis = findgen((maxr23-minr23)/dr23+1)*dr23+minr23
    maxo32 = 0.5 & mino32 = -1.0 & do32 = 0.01
    o32axis = findgen((maxo32-mino32)/do32+1)*do32+mino32

    x = r23axis & y = o32axis

    biglogoh_upper = 9.72D - 0.777D*x - 0.951D*x^2 - 0.072D*x^3 - 0.811D*x^4 - $
      logq[0]*(0.0737D - 0.0713D*x - 0.141D*x^2 + 0.0373D*x^3 - 0.058D*x^4)
    biglogoh_lower = 9.40D + 4.65D*x - 3.17D*x^2 - logq[0]*(0.272D + 0.547D*x - 0.513D*x^2)

    keep = where(biglogoh_upper ge biglogoh_lower)
    biglogoh_upper = biglogoh_upper[keep]
    biglogoh_lower = biglogoh_lower[keep]
    r23axis = r23axis[keep]

    if (logoh lt max(biglogoh_lower)) then biglogoh = biglogoh_lower else biglogoh = biglogoh_upper
    r23 = interpol(r23axis,biglogoh,logoh)
    
    biglogq = (32.81D - 1.153*y^2 + logoh*(-3.396 - 0.025*y + 0.1444*y^2)) / $
      (4.603D - 0.3119*y - 0.163*y^2 + logoh*(-0.48 + 0.0271*y + 0.02037*y^2))
    o32 = interpol(o32axis,biglogq,logq)

    if keyword_set(debug) then begin
       djs_plot, r23axis, biglogoh_upper, ps=4, yrange=[7.0,9.3], xsty=3, ysty=3, $
         charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, xtitle='log R_{23}', $
         ytitle='12 + log (O/H)', title='12+log(O/H) = '+string(logoh,format='(F4.2)')+$
         ', log q = '+string(logq,format='(F4.2)')
       oplot, r23axis, biglogoh_lower, ps=4
       oplot, !x.crange, logoh*[1,1], line=0
       oplot, r23*[1,1], !y.crange, line=0
       cc = get_kbrd(1)
    endif

return, r23
end
    
